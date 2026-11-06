#!/usr/bin/env python3
"""Oprava CRS v COG bez reprojekcie (Python 3.9+, GDAL 3.8+ a osgeo_utils).

Pouzitie: python fix_cog_crs.py [PRIECINOK_ALEBO_TIFF] --crs EPSG:8353 [--check]
          python fix_cog_crs.py [PRIECINOK_ALEBO_TIFF] --reference vzor.tif
          python fix_cog_crs.py [PRIECINOK_ALEBO_TIFF] --crs EPSG:8353 --skip-data-check
          python fix_cog_crs.py [PRIECINOK_ALEBO_TIFF] --crs EPSG:8353 --metadata-only

Predvoleny priecinok je vedla skriptu; spracuju sa .tif/.tiff bez podpriecinkov.
CRS musi urcit pouzivatel: z parametrov projekcie sa datum neda spolahlivo odhadnut.
--reference prebera CRS referencneho rastra, nie jeho mriezku ani pixely.
Pouzite jednotky a suradnice uz musia zodpovedat zadanemu CRS.

V predvolenom rezime sa original neotvara na zapis: CRS sa priradi kopii VRT a z nej sa
vytvori novy COG cez GDAL COG CreateCopy. Vystup je BigTIFF s bezstratovou
kompresiou DEFLATE; velkost a kompresia sa mozu oproti originalu zmenit.
Existujuce prehlady sa kopiruju bez prepocitania; nove sa negeneruju.
Pred nahradenim sa overia presne pixely pasiem, masiek a prehladov,
mriezka, NoData, metadata, CRS a struktura COG (plna validacia osgeo_utils).
--skip-data-check vynecha hashe pixelov pasiem, masiek aj prehladov a plnu
kontrolu datovych blokov COG. Zachova kontrolu metadat a zakladnej struktury COG,
ale neoveruje zhodu hodnot rastra. Vytvorenie a kompresia celeho COG zostavaju.
--metadata-only je este rychlejsi: na mieste prepise iba GeoKeys v existujucom
priestore CRS poli. IFD, velkost suboru, komprimovane obrazove bloky a ich posuny
sa nemenia. Bez kopirovania, kompresie alebo hashovania rastra. Ak sa nova definicia
nezmesti alebo su metadata nejednoznacne, oprava sa odmietne pred zapisom.
--metadata-only --check overi vhodnost bez zapisu. Pri chybe naslednej kontroly
sa obnovia povodne bajty CRS z pamate. Pri pade procesu/vypadku napajania pocas
priameho zapisu nie je obnova zarucena; tento rezim nema atomicke nahradenie.
Docasny vystup je vedla originalu: treba miesto na novy COG.
Trvala zaloha sa nevytvara. Vstup musi byt COG, bez externych .aux.xml/.ovr/.msk
a bez GCP/viacerych obrazov. Subor nesmie byt otvoreny v QGIS/inom programe.
Tento skript je samostatny; fix_tif_crs.py nepotrebuje.
"""

import argparse
import hashlib
import math
import os
from pathlib import Path
import shutil
import stat
import struct
import sys
import tempfile
import io
import uuid

try:
    from osgeo import gdal, osr
except ImportError:
    gdal = osr = None

try:
    from osgeo_utils.samples.validate_cloud_optimized_geotiff import validate as validate_cog
except ImportError:
    validate_cog = None


def open_tif(path):
    dataset = gdal.OpenEx(str(path), gdal.OF_RASTER)
    if dataset is None:
        raise ValueError(f"Raster sa nepodarilo otvorit: {path}")
    if dataset.GetDriver().ShortName != "GTiff":
        dataset.Close()
        raise ValueError("Podporovany je iba GeoTIFF")
    return dataset


def target_crs(crs, reference):
    if reference:
        with gdal.OpenEx(str(reference), gdal.OF_RASTER) as dataset:
            spatial_ref = dataset.GetSpatialRef()
            if spatial_ref is None:
                raise ValueError("Referencny raster nema CRS")
            result = spatial_ref.Clone()
    else:
        result = osr.SpatialReference()
        if result.SetFromUserInput(crs) != 0:
            raise ValueError("Neplatny cielovy CRS")
    if result.Validate() != 0:
        raise ValueError("Cielovy CRS nepresiel validaciou GDAL")
    return result


def same_crs(left, right):
    # TIFF pouziva tradicne GIS poradie suradnic; samotny CRS ostava nezmeneny.
    if left is None or right is None:
        return False
    left, right = left.Clone(), right.Clone()
    left.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    right.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return bool(left.IsSame(right))


def pixel_hash(band):
    digest = hashlib.sha256()
    # Kanonicke poradie po riadkoch nezavisi od fyzickej velkosti TIFF blokov.
    # Bez konverzie datoveho typu, vratane presneho zapisu float/NaN.
    bytes_per_pixel = max(1, (gdal.GetDataTypeSize(band.DataType) + 7) // 8)
    rows = max(1, (8 * 1024 * 1024) // (band.XSize * bytes_per_pixel))
    for y in range(0, band.YSize, rows):
        data = band.ReadRaster(0, y, band.XSize, min(rows, band.YSize - y))
        if data is None:
            raise ValueError("Nepodarilo sa precitat pixely")
        digest.update(data)
    return digest.hexdigest()


def metadata(obj):
    result = {domain: obj.GetMetadata(domain) for domain in obj.GetMetadataDomainList() or []
              if domain not in ("DERIVED_SUBDATASETS", "IMAGE_STRUCTURE")}
    if isinstance(obj, gdal.Dataset):
        # Chybajuci RasterType GeoKey znamena PixelIsArea. GDAL po zapise CRS
        # moze tento predvoleny vyznam explicitne vystavit v metadatach.
        result.setdefault("", {}).setdefault("AREA_OR_POINT", "Area")
    return result


def snapshot(dataset, check_pixels=True):
    """Ne-CRS metadata a volitelne presne hashe pixelov, masiek a prehladov."""
    hash_pixels = pixel_hash if check_pixels else lambda band: None
    bands = []
    for i in range(1, dataset.RasterCount + 1):
        band = dataset.GetRasterBand(i)
        nodata = band.GetNoDataValue()
        if isinstance(nodata, float) and math.isnan(nodata):
            nodata = "NaN"
        table = band.GetColorTable()
        bands.append((band.DataType, band.GetDescription(), nodata,
                      band.GetScale(), band.GetOffset(), band.GetUnitType(),
                      band.GetColorInterpretation(), band.GetCategoryNames(),
                      [table.GetColorEntry(j) for j in range(table.GetCount())] if table else None,
                      metadata(band), hash_pixels(band), band.GetMaskFlags(),
                      hash_pixels(band.GetMaskBand()),
                      [(band.GetOverview(j).XSize, band.GetOverview(j).YSize,
                        hash_pixels(band.GetOverview(j)),
                        hash_pixels(band.GetOverview(j).GetMaskBand()))
                       for j in range(band.GetOverviewCount())]))
    return (dataset.RasterXSize, dataset.RasterYSize, dataset.RasterCount,
            dataset.GetGeoTransform(can_return_null=True), metadata(dataset), bands)


def signature(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ino


def tiff_layout(file):
    """Precita IFD a metadatove posuny; necita obrazove bloky."""
    size = file.seek(0, 2)

    def read(offset, length):
        if offset < 0 or length < 0 or offset + length > size or length > 64 * 1024 * 1024:
            raise ValueError("Neplatny alebo prilis velky metadatovy blok TIFF")
        file.seek(offset)
        data = file.read(length)
        if len(data) != length:
            raise ValueError("Neuplne metadata TIFF")
        return data

    header = read(0, 8)
    if header[:2] not in (b"II", b"MM"):
        raise ValueError("Neplatna hlavicka TIFF")
    endian = "<" if header[:2] == b"II" else ">"
    magic = struct.unpack(endian + "H", header[2:4])[0]
    if magic == 42:
        count_fmt, number_fmt, entry_fmt, width, header_size = "H", "I", "HHII", 4, 8
        offset = struct.unpack(endian + "I", header[4:8])[0]
    elif magic == 43 and struct.unpack(endian + "HH", header[4:8]) == (8, 0):
        count_fmt, number_fmt, entry_fmt, width, header_size = "Q", "Q", "HHQQ", 8, 16
        offset = struct.unpack(endian + "Q", read(8, 8))[0]
    else:
        raise ValueError("Nepodporovany format TIFF")
    type_sizes = {1: 1, 2: 1, 3: 2, 4: 4, 5: 8, 6: 1, 7: 1, 8: 2,
                  9: 4, 10: 8, 11: 4, 12: 8, 13: 4, 16: 8, 17: 8, 18: 8}
    count_size = struct.calcsize(count_fmt)
    entry_size = struct.calcsize(endian + entry_fmt)
    directories, protected, seen = [], [(0, header_size)], set()
    while offset:
        if offset in seen or len(seen) >= 4096:
            raise ValueError("Cyklicke alebo prilis pocetne IFD")
        seen.add(offset)
        count = struct.unpack(endian + count_fmt, read(offset, count_size))[0]
        if count > 4096:
            raise ValueError("Prilis vela poloziek IFD")
        end = offset + count_size + count * entry_size + width
        protected.append((offset, end))
        tags = {}
        for i in range(count):
            start = offset + count_size + i * entry_size
            tag, kind, elements, value = struct.unpack(endian + entry_fmt, read(start, entry_size))
            if tag in tags or kind not in type_sizes:
                raise ValueError("Nejednoznacna polozka TIFF")
            length = elements * type_sizes[kind]
            location = start + entry_size - width if length <= width else value
            if location + length > size:
                raise ValueError("Polozka TIFF presahuje subor")
            tags[tag] = (kind, elements, location, length, length > width)
        if 330 in tags:
            raise ValueError("Rezim --metadata-only nepodporuje SubIFD")
        # Chranime aj komprimovane obrazove data vsetkych prehladov a masiek.
        for offsets_tag, lengths_tag in ((273, 279), (324, 325)):
            if offsets_tag not in tags and lengths_tag not in tags:
                continue
            if offsets_tag not in tags or lengths_tag not in tags:
                raise ValueError("Chybaju posuny alebo dlzky obrazovych blokov")
            arrays = []
            for tag in (offsets_tag, lengths_tag):
                kind, elements, location, length, _ = tags[tag]
                if kind not in (3, 4, 16):
                    raise ValueError("Neplatny typ posunov obrazovych blokov")
                arrays.append(struct.unpack(endian + str(elements) + {3: "H", 4: "I", 16: "Q"}[kind],
                                            read(location, length)))
            if len(arrays[0]) != len(arrays[1]):
                raise ValueError("Nesuhlasny pocet obrazovych blokov")
            for start, length in zip(*arrays):
                if start + length > size:
                    raise ValueError("Obrazovy blok presahuje subor")
                if length:
                    protected.append((start, start + length))
        directories.append(tags)
        offset = struct.unpack(endian + number_fmt, read(end - width, width))[0]
    if not directories:
        raise ValueError("TIFF nema IFD")
    return endian, directories, protected, read


def crs_patches(path, desired, area_or_point):
    """GDAL vygeneruje GeoKeys; planuje sa iba prepis existujucich CRS poli."""
    with path.open("rb") as file:
        endian, directories, protected, read = tiff_layout(file)
        if any(any(tag in tags for tag in (34735, 34736, 34737)) for tags in directories[1:]):
            raise ValueError("Rezim --metadata-only nepodporuje CRS vo viacerych IFD")
        template = "/vsimem/fix_crs_" + uuid.uuid4().hex + ".tif"
        try:
            with gdal.GetDriverByName("GTiff").Create(template, 1, 1, 1, options=[
                    "ENDIANNESS=" + ("LITTLE" if endian == "<" else "BIG")]) as ds:
                ds.SetGeoTransform((0, 1, 0, 0, 0, -1))
                ds.SetMetadataItem("AREA_OR_POINT", area_or_point or "Area")
                if ds.SetSpatialRef(desired) != 0:
                    raise ValueError("Generovanie GeoKeys zlyhalo")
            memory = io.BytesIO(bytes(gdal.VSIGetMemFileBuffer_unsafe(template)))
        finally:
            gdal.Unlink(template)
        _, target_ifds, _, target_read = tiff_layout(memory)
        generated_values = {tag: target_read(value[2], value[3]) for tag, value in target_ifds[0].items()
                            if tag in (34735, 34736, 34737)}
        directory = generated_values.get(34735, b"")
        old_directory = directories[0].get(34735)
        if old_directory and len(directory) > old_directory[3]:
            words = struct.unpack(endian + str(len(directory) // 2) + "H", directory)
            entries = [words[i:i + 4] for i in range(4, 4 + words[3] * 4, 4)]
            # CitationGeoKeys su iba volitelne popisy. EPSG a vsetky definicne
            # parametre zostavaju. Samoodkazove SHORT hodnoty by potrebovali
            # preindexovanie, preto sa takyto variant nekompaktuje.
            if len(words) == 4 + words[3] * 4 and not any(key[1] == 34735 for key in entries):
                entries = [key for key in entries if key[0] not in (1026, 2049, 3073, 4097)]
                # Standardny dvojrozmerny EPSG CRS uz urcuje datum, projekciu
                # aj jednotky. OGC GeoTIFF 1.1, priklady F.2.1/F.2.2, povoluje
                # tento minimalny zapis. Nikdy ho nepouzivame pre vlastny CRS.
                authority = desired.GetAuthorityCode(None)
                if desired.GetAuthorityName(None) == "EPSG" and authority and desired.GetAxesCount() == 2:
                    crs_key = 3072 if desired.IsProjected() else 2048 if desired.IsGeographic() else None
                    matching = any(key == (crs_key, 0, 1, int(authority)) for key in entries)
                    if matching and 0 < int(authority) < 32767:
                        entries = [key for key in entries if key[0] in (1024, 1025, crs_key)]
                compact = words[:3] + (len(entries),) + tuple(value for key in entries for value in key)
                generated_values[34735] = struct.pack(endian + str(len(compact)) + "H", *compact)
                if not any(key[1] == 34737 for key in entries):
                    generated_values.pop(34737, None)
        patches = []
        for tag in (34735, 34736, 34737):
            if tag not in generated_values:
                continue  # Nepouzite stare GeoDouble/Ascii hodnoty uz nove GeoKeys nereferencuju.
            new_kind = target_ifds[0][tag][0]
            generated = generated_values[tag]
            new_length = len(generated)
            old = directories[0].get(tag)
            if old is None or old[0] != new_kind or not old[4] or new_length > old[3]:
                raise ValueError(f"Novy CRS sa nezmesti do existujucich metadat (tag {tag}, "
                                 f"potrebne {new_length} B, dostupne {old[3] if old else 0} B); "
                                 "pouzite rezim bez --metadata-only")
            _, _, start, length, _ = old
            # Nijaky zapis nesmie zasiahnut IFD, ine pole ani obrazove data.
            other_fields = [(value[2], value[2] + value[3]) for tags in directories
                            for key, value in tags.items() if value[4] and not (tags is directories[0] and key == tag)]
            if any(start < right and left < start + length for left, right in protected + other_fields):
                raise ValueError("CRS pole sa prekryva s inymi datami TIFF")
            original = read(start, length)
            if new_kind == 2:
                # TIFF ASCII smie mat NUL az na konci celeho pola. GeoKeys
                # referencuju iba skutocne retazce, nie vyplnovacie medzery.
                text = generated.rstrip(b"\0")
                corrected = text + b" " * (length - len(text) - 1) + b"\0"
            else:
                corrected = generated + b"\0" * (length - new_length)
            if original != corrected:
                patches.append((start, original, corrected))
        if 34735 not in target_ifds[0] or 34735 not in directories[0]:
            raise ValueError("Chyba GeoKeyDirectoryTag")
        return patches


def patch_crs_only(path, desired, original_signature, before, patches, check):
    if check:
        print(f"  NA OPRAVU: {sum(len(p[1]) for p in patches)} bajtov metadat; bez zapisu", flush=True)
        return True
    with path.open("r+b") as file:
        current = os.fstat(file.fileno())
        if (current.st_size, current.st_mtime_ns, current.st_ino) != original_signature:
            raise ValueError("Subor sa pocas pripravy zmenil")
        for offset, old, _ in patches:
            file.seek(offset)
            if file.read(len(old)) != old:
                raise ValueError("CRS metadata sa pocas pripravy zmenili")
        try:
            for offset, _, new in patches:
                file.seek(offset)
                if file.write(new) != len(new):
                    raise OSError("Neuplny zapis CRS")
            file.flush()
            os.fsync(file.fileno())
            with open_tif(path) as dataset:
                if not same_crs(dataset.GetSpatialRef(), desired):
                    raise ValueError("Zapisany CRS nesuhlasi")
                if snapshot(dataset, check_pixels=False) != before:
                    raise ValueError("Zmenili sa ne-CRS metadata")
                if dataset.GetMetadataItem("LAYOUT", "IMAGE_STRUCTURE") != "COG":
                    raise ValueError("Zmenilo sa usporiadanie COG")
                _, errors, _ = validate_cog(dataset, full_check=False)
                if errors:
                    raise ValueError("Validacia COG zlyhala: " + "; ".join(errors))
            if os.fstat(file.fileno()).st_size != current.st_size:
                raise ValueError("Zmenila sa velkost suboru")
            for offset, _, new in patches:
                file.seek(offset)
                if file.read(len(new)) != new:
                    raise ValueError("Kontrola zapisu CRS zlyhala")
        except BaseException:
            for offset, old, _ in patches:
                file.seek(offset)
                file.write(old)
            file.flush()
            os.fsync(file.fileno())
            raise
    print("  OPRAVENE: iba CRS metadata; obrazove bloky a ich posuny neboli zapisovane", flush=True)
    return True


def process(path, desired, check=False, skip_data_check=False, metadata_only=False):
    if path.is_symlink():
        raise ValueError("Symbolicke odkazy sa nespracuvaju")
    for suffix in (".aux.xml", ".ovr", ".msk"):
        if Path(str(path) + suffix).exists():
            raise ValueError(f"TIFF ma externy subor {suffix}; tento variant nie je podporovany")
    original = signature(path)
    with open_tif(path) as dataset:
        if dataset.GetSubDatasets() or dataset.GetMetadata("SUBDATASETS"):
            raise ValueError("Viacobrazkovy TIFF nie je podporovany")
        if dataset.GetGCPCount():
            raise ValueError("Raster s GCP nie je podporovany")
        if dataset.GetMetadataItem("LAYOUT", "IMAGE_STRUCTURE") != "COG":
            raise ValueError("Vstup nie je rozpoznany ako COG; pre obycajny TIFF pouzite fix_tif_crs.py")
        block_x, block_y = dataset.GetRasterBand(1).GetBlockSize()
        if block_x != block_y or block_x % 16:
            raise ValueError("Nepodporovana velkost blokov vstupneho COG")
        current = dataset.GetSpatialRef()
        print(f"{path.name}: {current.GetName() if current else 'BEZ CRS'} -> {desired.GetName()}", flush=True)
        if same_crs(current, desired):
            return False
        if metadata_only:
            before = snapshot(dataset, check_pixels=False)
            area_or_point = dataset.GetMetadataItem("AREA_OR_POINT")
        else:
            before = None
        if not metadata_only:
            if check:
                print("  NA OPRAVU (bez zapisu)", flush=True)
                return True
            if skip_data_check:
                print("  VYTVORENIE NOVEHO COG (bez kontroly zhody pixelov)", flush=True)
            else:
                print("  KONTROLA PIXELOV A VYTVORENIE NOVEHO COG", flush=True)
            before = snapshot(dataset, check_pixels=not skip_data_check)

    if metadata_only:
        patches = crs_patches(path, desired, area_or_point)
        return patch_crs_only(path, desired, original, before, patches, check)

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(prefix=".fix_cog_crs_", suffix=".tif",
                                         dir=path.parent, delete=False) as file:
            temporary = Path(file.name)
        with open_tif(path) as dataset:
            # VRT je iba v pamati; odkazuje na povodne pixely aj ich prehlady.
            # Priradenie CRS nemení zdroj ani nepouziva reprojekciu.
            with gdal.GetDriverByName("VRT").CreateCopy("", dataset) as virtual:
                if virtual.SetSpatialRef(desired) != 0:
                    raise ValueError("Priradenie CRS zlyhalo")
                options = ["COMPRESS=DEFLATE", "BIGTIFF=YES", f"BLOCKSIZE={block_x}",
                           "OVERVIEWS=FORCE_USE_EXISTING", "OVERVIEW_COMPRESS=DEFLATE"]
                with gdal.GetDriverByName("COG").CreateCopy(str(temporary), virtual,
                                                            strict=1, options=options):
                    pass
        with open_tif(temporary) as dataset:
            if not same_crs(dataset.GetSpatialRef(), desired):
                raise ValueError("Zapisany CRS nesuhlasi s pozadovanym")
            if snapshot(dataset, check_pixels=not skip_data_check) != before:
                if skip_data_check:
                    raise ValueError("Zmenila sa mriezka, parametre pasiem/prehladov alebo ne-CRS metadata")
                raise ValueError("Zmenili sa pixely, masky, prehlady alebo ne-CRS metadata")
            if dataset.GetMetadataItem("LAYOUT", "IMAGE_STRUCTURE") != "COG":
                raise ValueError("Vystup nema usporiadanie COG")
            _, errors, _ = validate_cog(dataset, full_check=not skip_data_check)
            if errors:
                raise ValueError("Vystup nepresiel validaciou COG: " + "; ".join(errors))
        if Path(str(temporary) + ".aux.xml").exists():
            raise ValueError("GDAL ulozil metadata externe namiesto do TIFF")
        if signature(path) != original:
            raise ValueError("Original sa pocas opravy zmenil")
        with temporary.open("r+b") as file:
            os.fsync(file.fileno())
        shutil.copymode(path, temporary)
        try:
            os.replace(temporary, path)
        except PermissionError as error:
            raise PermissionError("Nahradenie COG odmietnute; original zostal nezmeneny. "
                                  "Zavrite subor v QGIS/inom programe a skontrolujte "
                                  "opravnenia a atribut iba na citanie.") from error
        if skip_data_check:
            print("  OPRAVENE; CRS, metadata a zakladna struktura COG overene; zhoda pixelov nekontrolovana", flush=True)
        else:
            print("  OPRAVENE; COG overeny, pixely, masky, prehlady a mriezka nezmenene", flush=True)
        return True
    finally:
        for temp in (temporary,):
            if temp is not None:
                for item in (temp, Path(str(temp) + ".aux.xml")):
                    if item.exists():
                        item.chmod(item.stat().st_mode | stat.S_IWUSR)
                        item.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--crs", help="Spravny CRS, napr. EPSG:8353 (neodhaduje sa automaticky)")
    source.add_argument("--reference", type=Path, help="Raster so spravnym CRS")
    parser.add_argument("--check", action="store_true", help="Iba kontrola, bez zapisu a hashovania pixelov")
    parser.add_argument("--skip-data-check", action="store_true",
                        help="Vynechat kontrolne citanie/hashovanie pixelov a plnu validaciu blokov COG; "
                             "kontrola CRS, mriezky, metadat a zakladnej struktury zostava")
    parser.add_argument("--metadata-only", action="store_true",
                        help="Najrychlejsie: prepisat iba existujuce CRS polia na mieste, bez kopirovania "
                             "rastra. Ak sa novy CRS nezmesti, subor sa neupravi. --check iba overi vhodnost.")
    args = parser.parse_args(argv)
    if gdal is None:
        parser.error("Chyba GDAL Python (osgeo). Spustite napr. cez QGIS python-qgis-ltr.bat.")
    if validate_cog is None:
        parser.error("Chyba osgeo_utils.samples.validate_cloud_optimized_geotiff z balika GDAL.")
    if int(gdal.VersionInfo("VERSION_NUM")) < 3080000 or gdal.GetDriverByName("COG") is None:
        parser.error("Vyžaduje GDAL 3.8+ s ovladacom COG.")
    gdal.UseExceptions()
    osr.UseExceptions()
    try:
        desired = target_crs(args.crs, args.reference)
    except (RuntimeError, ValueError, OSError) as error:
        parser.error(str(error))
    if args.path.is_dir():
        files = sorted(p for p in args.path.iterdir() if p.is_file() and p.suffix.lower() in (".tif", ".tiff"))
    elif args.path.is_file() and args.path.suffix.lower() in (".tif", ".tiff"):
        files = [args.path]
    else:
        parser.error("Zadajte existujuci priecinok alebo TIFF")
    changed = errors = 0
    for path in files:
        try:
            changed += process(path, desired, args.check, args.skip_data_check, args.metadata_only)
        except (OSError, ValueError, RuntimeError, struct.error) as error:
            errors += 1
            print(f"CHYBA: {path.name}: {error}", file=sys.stderr, flush=True)
    label = "na opravu" if args.check else "opravene"
    print(f"Spolu: {len(files)}, {label}: {changed}, bez zmeny: {len(files)-changed-errors}, chyby: {errors}.", flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nPrerusene pouzivatelom.", file=sys.stderr)
        sys.exit(130)
