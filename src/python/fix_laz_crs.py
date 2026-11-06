#!/usr/bin/env python3
"""Doplnenie chybajuceho WKT do LAZ bez reprojekcie (Python 3.9+).

Pouzitie: python fix_laz_crs.py [PRIECINOK_ALEBO_LAZ] --crs EPSG:8353 [--check]
Instalacia: python -m pip install "laspy[lazrs,pyproj]>=2.5,<3"

Predvolene spracuje LAZ vedla skriptu, bez podpriecinkov. Neprázdny WKT
vo VLR/EVLR (LASF_Projection aj liblas) preskoci. Existujuce GeoKeys musia
zodpovedat zadanemu CRS; nezhoda alebo necitatelny CRS je chyba.
Suradnice uz musia byt v zadanom CRS. Skript ich netransformuje.

Zapis prebieha cez docasny LAZ v rovnakom priecinku a rekompresiu bodov.
Pred nahradenim originalu overi hash vsetkych dekomprimovanych bodov,
ich pocet, mierky, posuny a ostatne VLR/EVLR. Starsi LAS povysi na 1.4,
aby mohol zapisat standardny WKT; format bodov zostane zachovany.
Vyzaduje volne miesto na dalsiu kopiu LAZ, trvalu zalohu nevytvara.
COPC a waveform odmieta. Subory pocas behu nesmu byt zapisovane inym programom.
"""

import argparse
import copy
import hashlib
import os
from pathlib import Path
import sys
import tempfile

try:
    import laspy
    from laspy.vlrs.vlrlist import VLRList
    from pyproj import CRS
except ImportError:
    laspy = CRS = None


PROJECTION_USERS = {"LASF_Projection", "liblas"}
CRS_RECORDS = {2111, 2112, 34735, 34736, 34737}
CHUNK_POINTS = 500_000


def is_projection(record, ids=CRS_RECORDS):
    return record.user_id in PROJECTION_USERS and record.record_id in ids


def records(header):
    return list(header.vlrs) + list(header.evlrs or [])


def has_wkt(header):
    return any(is_projection(v, {2112}) and v.record_data_bytes().rstrip(b"\0").strip()
               for v in records(header))


def check_input(header, crs):
    """Nikdy neprepisuje iny znamy CRS ani neignoruje necitatelne GeoKeys."""
    if has_wkt(header):
        return False
    if any(v.user_id.lower() == "copc" for v in records(header)):
        raise ValueError("COPC nie je podporovany")
    if header.point_format.id in (4, 5, 9, 10) or header.start_of_waveform_data_packet_record:
        raise ValueError("Waveform data nie su podporovane")
    for v in records(header):
        if is_projection(v, {34735, 34736, 34737}):
            # liblas GeoKeys nie su standardny vstup parse_crs: neodhadovat ich vyznam.
            if v.user_id != "LASF_Projection":
                raise ValueError("Nestandardne GeoKeys liblas; vyzaduju manualnu kontrolu")
    if any(is_projection(v, {34735, 34736, 34737}) for v in records(header)):
        known = header.parse_crs(prefer_wkt=False)
        if known is None:
            raise ValueError("Existujuce GeoKeys sa nepodarilo precitat")
        if not known.equals(crs, ignore_axis_order=True):
            raise ValueError("Existujuce GeoKeys obsahuju iny CRS ako --crs")
    return True


def signature(path):
    stat = path.stat()
    return stat.st_size, stat.st_mtime_ns, stat.st_ino


def metadata(header):
    def other(vlrs):
        return [(v.user_id, v.record_id, v.description, v.record_data_bytes())
                for v in vlrs or [] if not is_projection(v)
                and not (v.user_id == "laszip encoded" and v.record_id == 22204)]
    # Kompresny VLR regeneruje backend. Ostatne zaznamy musia ostat identicke.
    return (header.point_format.id, header.point_format.size,
            tuple(header.scales), tuple(header.offsets), header.point_count,
            tuple(header.mins), tuple(header.maxs), header.file_source_id,
            header.uuid, header.system_identifier, header.generating_software,
            header.creation_date, header.global_encoding.value & ~16,
            header.extra_header_bytes, header.extra_vlr_bytes,
            other(header.vlrs), other(header.evlrs))


def verify(path, expected_hash, expected_metadata, crs):
    with laspy.open(path) as reader:
        actual_metadata = metadata(reader.header)
        if actual_metadata != expected_metadata:
            raise ValueError("Overenie ostatnych metadat zlyhalo")
        actual_crs = reader.header.parse_crs()
        if actual_crs is None or not actual_crs.equals(crs):
            raise ValueError("Overenie zapisaneho CRS zlyhalo")
        if not has_wkt(reader.header) or not reader.header.global_encoding.wkt:
            raise ValueError("Vystup neobsahuje oznaceny WKT")
        digest = hashlib.sha256()
        count = 0
        for points in reader.chunk_iterator(CHUNK_POINTS):
            digest.update(points.array.tobytes())
            count += len(points)
        if digest.digest() != expected_hash or count != reader.header.point_count:
            raise ValueError("Overenie vsetkych bodov zlyhalo")


def repair(path, crs, check=False):
    initial = signature(path)
    temporary = None
    try:
        with laspy.open(path) as reader:
            if not check_input(reader.header, crs):
                return False
            if check:
                return True
            original_metadata = metadata(reader.header)
            header = copy.deepcopy(reader.header)
            header.set_version_and_point_format(laspy.header.Version(1, 4), header.point_format)
            # Odstranit aj prazdne/duplicitne CRS zaznamy v EVLR a liblas.
            # Surove VLR zabránia laspy prepocitat aj minima/maxima Extra Bytes.
            # Zoznam menit na mieste: setter vlrs regeneruje Extra Bytes VLR.
            header.vlrs[:] = [laspy.VLR(user_id=v.user_id, record_id=v.record_id,
                                       description=v.description, record_data=v.record_data_bytes())
                              for v in reader.header.vlrs if not is_projection(v)
                              and not (v.user_id == "laszip encoded" and v.record_id == 22204)]
            header.evlrs = VLRList(v for v in header.evlrs or [] if not is_projection(v))
            header.add_crs(crs, keep_compatibility=False)
            fd, name = tempfile.mkstemp(prefix=".fix_laz_crs_", suffix=".laz", dir=path.parent)
            os.close(fd)
            temporary = Path(name)
            digest = hashlib.sha256()
            count = 0
            with laspy.open(temporary, mode="w", header=header, do_compress=True) as writer:
                for points in reader.chunk_iterator(CHUNK_POINTS):
                    digest.update(points.array.tobytes())
                    count += len(points)
                    writer.write_points(points)
                if header.evlrs:
                    writer.write_evlrs(header.evlrs)
            if count != reader.header.point_count:
                raise ValueError("Pocet precitanych bodov nesuhlasi s hlavickou")
        verify(temporary, digest.digest(), original_metadata, crs)
        if signature(path) != initial:
            raise ValueError("Povodny subor sa pocas spracovania zmenil")
        with temporary.open("r+b") as file:
            file.flush()
            os.fsync(file.fileno())
        os.replace(temporary, path)
        return True
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--crs", required=True, help="CRS bodov, napr. EPSG:8353")
    parser.add_argument("--check", action="store_true", help="Iba kontrola, bez zapisu a dekompresie bodov")
    args = parser.parse_args(argv)
    if laspy is None or CRS is None:
        parser.error('Chybaju kniznice. Instalacia: python -m pip install "laspy[lazrs,pyproj]>=2.5,<3"')
    try:
        crs = CRS.from_user_input(args.crs)
    except Exception as error:
        parser.error(f"Neplatny CRS: {error}")
    if args.path.is_symlink():
        parser.error("Symbolicke odkazy nie su podporovane")
    if args.path.is_dir():
        files = sorted(p for p in args.path.iterdir() if p.is_file()
                       and p.suffix.lower() == ".laz" and not p.is_symlink())
    elif args.path.is_file() and args.path.suffix.lower() == ".laz":
        files = [args.path]
    else:
        parser.error("Zadajte existujuci priecinok alebo LAZ subor")
    changed = skipped = errors = 0
    for path in files:
        try:
            print(f"KONTROLUJEM: {path.name}" if args.check else f"SPRACUVAM: {path.name}", flush=True)
            if repair(path, crs, check=args.check):
                changed += 1
                print(f"{'NA DOPLNENIE' if args.check else 'DOPLNENE A OVERENE'}: {path.name} ({args.crs})", flush=True)
            else:
                skipped += 1
                print(f"MA WKT, PRESKAKUJEM: {path.name}", flush=True)
        except Exception as error:
            errors += 1
            print(f"CHYBA: {path.name}: {error}", file=sys.stderr, flush=True)
    print(f"Spolu: {len(files)}, {'na doplnenie' if args.check else 'doplnene'}: {changed}, "
          f"preskocene: {skipped}, chyby: {errors}.", flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nPrerusene pouzivatelom.", file=sys.stderr)
        sys.exit(130)
