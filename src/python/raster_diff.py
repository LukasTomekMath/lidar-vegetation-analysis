import rasterio
from rasterio.windows import from_bounds
import numpy as np

# porovnanie dvoch rastrov na overlaped area, toto bolo v case study

r1 = r"vsetky_OLD\\Boky_OLD.tif"
r2 = r"vsetky_NEW\\Boky_NEW.tif"
out = "Boky_dif.tif"

OUT_NODATA = -9999.0

with rasterio.open(r1) as src1, rasterio.open(r2) as src2:
    print("OLD bands:", src1.count)
    print("NEW bands:", src2.count)

    n_bands = min(src1.count, src2.count)

    left = max(src1.bounds.left, src2.bounds.left)
    bottom = max(src1.bounds.bottom, src2.bounds.bottom)
    right = min(src1.bounds.right, src2.bounds.right)
    top = min(src1.bounds.top, src2.bounds.top)

    if left >= right or bottom >= top:
        raise ValueError("No overlapping area between rasters.")

    # OLD je ref grid
    win1 = from_bounds(left, bottom, right, top, src1.transform)
    win1 = win1.round_offsets().round_lengths()

    
    overlap_bounds = rasterio.windows.bounds(win1, src1.transform)

    win2 = from_bounds(*overlap_bounds, transform=src2.transform)
    win2 = win2.round_offsets().round_lengths()

    profile = src1.profile.copy()
    profile.update(
        count=n_bands,
        height=int(win1.height),
        width=int(win1.width),
        transform=src1.window_transform(win1),
        dtype="float32",
        nodata=OUT_NODATA
    )

    with rasterio.open(out, "w", **profile) as dst:
        for b in range(1, n_bands + 1):
            a = src1.read(
                b,
                window=win1,
                masked=True
            ).astype("float32")

            c = src2.read(
                b,
                window=win2,
                out_shape=a.shape,
                masked=True
            ).astype("float32")

            combined_mask = (
                np.ma.getmaskarray(a) |
                np.ma.getmaskarray(c)
            )

            diff = np.ma.array(c - a, mask=combined_mask)

            dst.write(
                diff.filled(OUT_NODATA).astype("float32"),
                b
            )

            desc = src1.descriptions[b - 1]
            if desc:
                dst.set_band_description(b, desc)

            valid = diff.compressed()
            if valid.size:
                dst.update_tags(
                    b,
                    STATISTICS_MINIMUM=float(valid.min()),
                    STATISTICS_MAXIMUM=float(valid.max()),
                    STATISTICS_MEAN=float(valid.mean()),
                    STATISTICS_STDDEV=float(valid.std())
                )

print(f"Done: {out}")