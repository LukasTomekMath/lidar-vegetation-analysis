import os
import re
import time
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from dotenv import load_dotenv
from datetime import date, timedelta
from pystac_client import Client
from pystac_client.stac_api_io import StacApiIO
from rasterio.errors import RasterioIOError, WindowError
from rasterio.windows import Window
from rasterio.mask import mask, raster_geometry_mask
from shapely.geometry import box, shape
from urllib3.util.retry import Retry

load_dotenv()

# ======================================================
DATE_START, DATE_END = "2024-04-01", "2024-04-14"

EXPAND_STEP_DAYS = 5                        # widen the date window by this much on each side per round (0 = never expand)
MAX_EXPAND_DAYS = 30                        # stop after this many days on each side

LIMIT = 2                                # None = all polygons
N_BEST = 1                                  # clearest scenes to extract per polygon
MIN_COVERAGE = 0.90                         # polygon must be >= 90% covered by valid tile data
MAX_CLOUD_PCT = 20                          # max % of INVALID_CLASSES in polygon
CATALOG_MAX_CLOUD_PCT = 20                  # max % of clouds in catalog
INVALID_CLASSES = [3, 7, 8, 9, 10, 11]      # cloud shadow, unclassified, medium/high prob cloud, thin cirrus, snow/ice
RETRIES = 5                                 # attempts per open/read before giving up

PAD_M = 40                                  # padding around bbox in meters
MAX_RESOLUTION = 20                         # coarsest band resolution in meters 

PREVIEW = False                             # plot the first extracted polygon at the end
CLIP_TO_BBOX = True                         # save the polygon's bounding box instead of the polygon shape
STATS_ON_BBOX = False                       # True = cloud % / coverage over the bbox instead of the polygon

OUT = Path("copernicus" + "_" + DATE_START + "_" + DATE_END)

BANDS_10M = ["B02_10m", "B03_10m", "B04_10m", "B08_10m"]
BANDS_20M = ["B05_20m", "B06_20m", "B07_20m", "B8A_20m", "B11_20m", "B12_20m", "SCL_20m"]
# ======================================================

OUT.mkdir(exist_ok=True)

access_key = os.getenv("CDSE_S3_ACCESS_KEY")
secret_key = os.getenv("CDSE_S3_SECRET_KEY")
SHP = os.getenv("SHP_FILEPATH")
if not access_key or not secret_key:
    raise RuntimeError("CDSE_S3_ACCESS_KEY / CDSE_S3_SECRET_KEY missing from .env")

os.environ.update(
    AWS_ACCESS_KEY_ID=access_key,
    AWS_SECRET_ACCESS_KEY=secret_key,
    AWS_S3_ENDPOINT="eodata.dataspace.copernicus.eu",
    AWS_VIRTUAL_HOSTING="FALSE",
    AWS_HTTPS="YES",
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    GDAL_HTTP_MAX_RETRY="10",      
    GDAL_HTTP_RETRY_DELAY="3",
    VSI_CACHE="TRUE",                
    VSI_CACHE_SIZE="67108864",
)

failures = []   # anything that still failed after retries

def safe(s):
    """Replace characters that are illegal in file names."""
    return re.sub(r'[<>:"/\\|?*]+', "-", str(s)).strip(" .") or "unnamed"


def vsis3(item, key):
    href = item.assets[key].href
    if not href.startswith("s3://"):
        raise ValueError(f"Expected an s3:// href, got: {href}")
    return "/vsis3/" + href.removeprefix("s3://")


def out_path(name, date, res):
    stem = name if N_BEST == 1 else f"{name}_{date}"   # date needed or scenes overwrite each other
    return OUT / f"{stem}_{res}.tif"


def with_retry(fn, tries=RETRIES, base_delay=3):
    """Retry transient S3 failures (429, 'Stream too short', ...) with growing delays."""
    for attempt in range(1, tries + 1):
        try:
            return fn()
        except ValueError:
            raise                       # e.g. polygon doesn't overlap the raster
        except Exception as e:
            if attempt == tries:
                raise
            wait = base_delay * attempt
            print(f"  retry {attempt}/{tries - 1} in {wait}s: {str(e)[:90]}")
            time.sleep(wait)


def open_src(path):
    src = rasterio.open(path)
    if src.transform.is_identity:       # header read was cut short: treat as transient
        src.close()
        raise RasterioIOError(f"No georeferencing read for {path}")
    return src


def open_all(item, keys):
    srcs = []
    try:
        for k in keys:
            srcs.append(open_src(vsis3(item, k)))
        return srcs
    except Exception:
        for s in srcs:
            s.close()
        raise


def read_scl(src, geom):
    outside, _, window = raster_geometry_mask(src, [geom], crop=True, all_touched=True)
    return outside, src.read(1, window=window)


def clip_all(srcs, geom):
    arrays = []
    for src in srcs:
        arr, tr = mask(src, [geom], crop=True, all_touched=True, nodata=0)
        arrays.append(arr[0])
    return np.stack(arrays), tr

def snapped_bounds(src, geom, pad_m = PAD_M, step = MAX_RESOLUTION):
    x0, y0 = src.transform.c, src.transform.f
    minx, miny, maxx, maxy = geom.bounds
    lo = lambda v, o: o + np.floor((v - o) / step) * step
    hi = lambda v, o: o + np.ceil((v - o) / step) * step
    return lo(minx - pad_m, x0), lo(miny - pad_m, y0), hi(maxx + pad_m, x0), hi(maxy + pad_m, y0)

def clip_bbox(srcs, geom):
    arrays = []
    for src in srcs:
        t = src.transform
        minx, miny, maxx, maxy = snapped_bounds(src, geom)
        win = Window(
            int(round((minx - t.c) / t.a)), int(round((t.f - maxy) / -t.e)),
            int(round((maxx - minx) / t.a)), int(round((maxy - miny) / -t.e)),
        )
        try:
            win = win.intersection(Window(0, 0, src.width, src.height))
        except WindowError:
            raise ValueError("no overlap with raster")   # handled as a skip, not retried
        arrays.append(src.read(1, window=win))
        tr = src.window_transform(win)
    return np.stack(arrays), tr


# load all polygons 
gdf = gpd.read_file(SHP, encoding="cp1250")
if LIMIT:
    gdf = gdf.head(LIMIT)
gdf["poly_id"] = gdf.index
gdf["name"] = (
    gdf["NAZOV"].fillna("").astype(str).str.strip() + "_" + gdf["KOD_PRALES"].fillna("").astype(str).str.strip()
).map(safe)
dup = gdf["name"].duplicated(keep=False)                 # same name twice -> add row index
gdf.loc[dup, "name"] = gdf.loc[dup, "name"] + "_" + gdf.loc[dup, "poly_id"].astype(str)
print(f"{len(gdf)} polygons ({dup.sum()} had duplicate names)")

# one search for the whole area, then match polygons to scenes 
gdf_4326 = gdf.to_crs(4326)
retry = Retry(
    total=8,
    backoff_factor=2,                        # growing waits between attempts
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=None,                    
    respect_retry_after_header=True,
)
stac = Client.open("https://stac.dataspace.copernicus.eu/v1/", stac_io=StacApiIO(max_retries=retry))

search_kwargs = {}
if CATALOG_MAX_CLOUD_PCT is not None:
    search_kwargs["query"] = {"eo:cloud_cover": {"lt": CATALOG_MAX_CLOUD_PCT}}


# ======================================================
# pass 1: SCL stats for every polygon-scene pair 
# ======================================================
def stats_for_item(item_id, pids):
    item = items_by_id[item_id]
    rows = []
    try:
        with rasterio.Env():
            src = with_retry(lambda: open_src(vsis3(item, "SCL_20m")))
            try:
                geoms = gdf.geometry.loc[pids].to_crs(src.crs)
                pixel_area = abs(src.res[0] * src.res[1])
                for pid, geom in zip(geoms.index, geoms):
                    geom = box(*snapped_bounds(src, geom)) if STATS_ON_BBOX else geom
                    try:
                        outside, scl = with_retry(lambda: read_scl(src, geom))
                    except ValueError:      # no overlap with this tile
                        continue
                    except Exception as e:
                        failures.append(dict(stage="pass1", poly_id=pid, item_id=item_id, res="20m", error=str(e)[:200]))
                        continue
                    valid = ~outside & (scl != 0)
                    n_valid = int(valid.sum())
                    rows.append(dict(
                        poly_id=pid,
                        item_id=item_id,
                        date=item.datetime.date(),
                        offset_days=(item.datetime.date() - center).days,
                        tile=item_id.split("_")[5],
                        scene_cloud=item.properties.get("eo:cloud_cover"),
                        valid_px=n_valid,
                        coverage=min(n_valid * pixel_area / geom.area, 1.0),
                        cloud_pct=100 * (np.isin(scl, INVALID_CLASSES) & valid).sum() / n_valid if n_valid else np.nan,
                    ))
            finally:
                src.close()
        print(f"pass 1: {item_id[:45]}... {len(rows)} polygons")
    except Exception as e:
        failures.append(dict(stage="pass1-open", poly_id=None, item_id=item_id, res="20m", error=str(e)[:200]))
        print(f"pass 1 FAILED for {item_id}: {str(e)[:100]}")
    return rows


d0, d1 = date.fromisoformat(DATE_START), date.fromisoformat(DATE_END)
center = d0 + timedelta(days=(d1 - d0).days // 2)

# round 0 = your window; later rounds = only the new strips on each side
rounds = [([(d0, d1)], 0)]                 # (date ranges, days of reach on each side)
reach = 0
while EXPAND_STEP_DAYS > 0 and reach < MAX_EXPAND_DAYS:
    nxt = min(reach + EXPAND_STEP_DAYS, MAX_EXPAND_DAYS)
    rounds.append((
        [(d0 - timedelta(days=nxt), d0 - timedelta(days=reach + 1)),
         (d1 + timedelta(days=reach + 1), d1 + timedelta(days=nxt))],
        nxt,
    ))
    reach = nxt

# ONE search covering the full widened range; the rounds below only filter it by date
max_reach = rounds[-1][1]
full_start, full_end = d0 - timedelta(days=max_reach), d1 + timedelta(days=max_reach)
all_items = with_retry(lambda: list(stac.search(
    collections=["sentinel-2-l2a"],
    bbox=list(gdf_4326.total_bounds),
    datetime=f"{full_start}/{full_end}",
    **search_kwargs,
).items()), base_delay=10)
if not all_items:
    raise SystemExit("The catalog search returned no scenes.")
items_by_id = {i.id: i for i in all_items}
all_scenes = gpd.GeoDataFrame(
    {"item_id": [i.id for i in all_items], "day": [i.datetime.date() for i in all_items]},
    geometry=[shape(i.geometry) for i in all_items],
    crs=4326,
)
print(f"{len(all_items)} scenes between {full_start} and {full_end}")

parts, n_pairs = [], 0
unresolved = set(gdf["poly_id"])
searched_days = 0

for n, (strips, reach) in enumerate(rounds):
    if not unresolved:
        break
    todo = gdf_4326[gdf_4326["poly_id"].isin(unresolved)]
    print(f"round {n}: {len(unresolved)} polygons unresolved, using scenes from " + " and ".join(f"{a} to {b}" for a, b in strips))

    in_round = np.zeros(len(all_scenes), dtype=bool)
    for a, b in strips:
        in_round |= ((all_scenes["day"] >= a) & (all_scenes["day"] <= b)).to_numpy()
    searched_days = reach
    new_scenes = all_scenes[in_round]
    if new_scenes.empty:
        continue

    pairs = gpd.sjoin(todo[["poly_id", "geometry"]], new_scenes[["item_id", "geometry"]], predicate="intersects")
    pairs = pd.DataFrame(pairs.drop(columns="geometry")).drop(columns="index_right", errors="ignore")
    n_pairs += len(pairs)
    by_item = {k: g["poly_id"].tolist() for k, g in pairs.groupby("item_id")}
    print(f"  {len(new_scenes)} scenes, {len(pairs)} polygon-scene pairs")

    results = [stats_for_item(i, by_item[i]) for i in by_item]
    rows = [r for rows_ in results for r in rows_]
    if rows:
        parts.append(pd.DataFrame(rows))
        so_far = pd.concat(parts, ignore_index=True)
        ok_now = (so_far["coverage"] >= MIN_COVERAGE) & (so_far["cloud_pct"] <= MAX_CLOUD_PCT)
        unresolved -= set(so_far.loc[ok_now, "poly_id"])

print(f"\n{len(unresolved)} polygons still unresolved after using scenes up to {searched_days} days around the window")

if not parts:
    raise SystemExit("No polygon-scene stats were computed. Check the errors above.")
index = pd.concat(parts, ignore_index=True)
print(f"pass 1: {len(index)} stats rows for {n_pairs} pairs, {len(failures)} failures")
index = index.merge(gdf[["poly_id", "name"]], on="poly_id")
index = index[["poly_id", "name"] + [c for c in index.columns if c not in ("poly_id", "name")]]
index.to_csv(OUT / "index_all.csv", index=False)

# choose the best scene(s) per polygon
good = index[(index["coverage"] >= MIN_COVERAGE) & (index["cloud_pct"] <= MAX_CLOUD_PCT)]
selected = (
    good.sort_values(["poly_id", "cloud_pct", "coverage"], ascending=[True, True, False])
    .groupby("poly_id")
    .head(N_BEST)
)
selected.to_csv(OUT / "selected.csv", index=False)

best_seen = index.groupby("poly_id").agg(best_cloud_pct=("cloud_pct", "min"), best_coverage=("coverage", "max"))
skipped = (
    gdf.loc[~gdf["poly_id"].isin(selected["poly_id"]), ["poly_id", "name"]]
    .merge(best_seen, on="poly_id", how="left")      # NaN = no scene overlapped at all
)
skipped.to_csv(OUT / "skipped.csv", index=False)
selected_num = selected['poly_id'].nunique()
print(f"\n{selected_num} polygons selected, {len(skipped)} skipped (see skipped.csv)")


# ======================================================
# pass 2: extract bands for the selected pairs 
# ======================================================
sel_by_item = {k: g for k, g in selected.groupby("item_id")}

total_written = 0
def extract_item(item_id):
    global total_written
    item = items_by_id[item_id]
    grp = sel_by_item[item_id]
    written = 0
    with rasterio.Env():
        for res, keys in (("10m", BANDS_10M), ("20m", BANDS_20M)):
            todo = [r for r in grp.itertuples() if not out_path(r.name, r.date, res).exists()]
            if not todo:
                continue
            try:
                srcs = with_retry(lambda: open_all(item, keys))
            except Exception as e:
                for r in todo:
                    failures.append(dict(stage="pass2-open", poly_id=r.poly_id, item_id=item_id, res=res, error=str(e)[:200]))
                print(f"pass 2: could not open {item_id[:45]}... ({res}): {str(e)[:90]}")
                continue
            try:
                crs = srcs[0].crs
                geoms = gdf.geometry.loc[[r.poly_id for r in todo]].to_crs(crs)
                for r in todo:
                    try:
                        clip = clip_bbox if CLIP_TO_BBOX else clip_all    # same (srcs, geom) signature
                        stack, tr = with_retry(lambda: clip(srcs, geoms.loc[r.poly_id]))
                    except ValueError as e:
                        print(f"  skip {r.name} ({res}): {e}")
                        continue
                    except Exception as e:
                        failures.append(dict(stage="pass2-read", poly_id=r.poly_id, item_id=item_id, res=res, error=str(e)[:200]))
                        print(f"  FAILED {r.name} ({res}): {str(e)[:90]}")
                        continue

                    path = out_path(r.name, r.date, res)
                    tmp = path.with_suffix(".tmp")      # write then rename, so a crash never leaves a half file
                    with rasterio.open(
                        tmp, "w", driver="GTiff",
                        height=stack.shape[1], width=stack.shape[2], count=stack.shape[0],
                        dtype=stack.dtype, crs=crs, transform=tr, compress="deflate", nodata=0,
                    ) as dst:
                        dst.write(stack)
                        dst.descriptions = tuple(k.rsplit("_", 1)[0] for k in keys)
                    os.replace(tmp, path)
                    total_written += 1
                    written += 1
                    print(f"pass 2 : wrote |{total_written} / {2*selected_num}| : {path.name}")
            finally:
                for s in srcs:
                    s.close()
    return written


total = sum(extract_item(item_id) for item_id in sel_by_item)
print(f"\nDone: {total} new files in {OUT}/")

if failures:
    pd.DataFrame(failures).to_csv(OUT / "failed.csv", index=False)
    print(f"{len(failures)} failures, see failed.csv. Re-run the script: existing files are skipped, "
          "so only the missing ones are retried.")
elif (OUT / "failed.csv").exists():
    (OUT / "failed.csv").unlink()


# ======================================================
# preview of the first extracted polygon
# ======================================================
if PREVIEW and len(selected):
    r = selected.iloc[0]
    p10, p20 = out_path(r["name"], r["date"], "10m"), out_path(r["name"], r["date"], "20m")
    if p10.exists() and p20.exists():
        with rasterio.open(p10) as s:
            rgb = s.read([3, 2, 1]).astype("float32")           # B04, B03, B02
        with rasterio.open(p20) as s:
            scl = s.read(s.count)                               # SCL is the last band
        rgb = np.clip((rgb - 1000) / 10000 / 0.3, 0, 1).transpose(1, 2, 0)
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        ax[0].imshow(rgb); ax[0].set_title(f"RGB (10 m) - {r['name']}")
        ax[1].imshow(scl, cmap="tab20", vmin=0, vmax=11); ax[1].set_title("SCL (20 m)")
        plt.show()