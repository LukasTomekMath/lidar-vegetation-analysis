#!/usr/bin/env python3
"""Oprava jednoduchych uvodzoviek vo WKT metadatách LAS/LAZ (Python 3.9+).

Bez argumentov spracuje *.laz vedla tohto skriptu, nie v pracovnom priecinku.
Pouzitie: python fix_laz_wkt_quotes.py [PRIECINOK] [--check]

Pouziva iba standardnu kniznicu Pythonu. Nevytvara zalohy ani protokoly.
Nemení CRS, dlzky zaznamov ani komprimovane body. Kontroluje syntax WKT,
nie geodeticku spravnost definicie. Chybajuce CRS nedopĺna.
Subory pocas behu nesmu byt zapisovane inym programom.
"""

import argparse
import hashlib
import os
from pathlib import Path
import re
import struct
import sys
from dataclasses import dataclass


MAX_WKT_BYTES = 16 * 1024 * 1024
CHUNK_BYTES = 8 * 1024 * 1024
WKT_ROOTS = {
    "PROJCS", "GEOGCS", "GEOCCS", "COMPD_CS", "VERT_CS", "LOCAL_CS",
    "PROJCRS", "GEOGCRS", "GEODCRS", "VERTCRS", "ENGCRS", "COMPOUNDCRS",
    "BOUNDCRS", "TIMECRS", "PARAMETRICCRS", "FITTED_CS",
}
TOKEN = re.compile(r'[A-Za-z_][A-Za-z_0-9]*|[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[\[\](),]')


@dataclass
class Patch:
    offset: int
    original: bytes
    corrected: bytes


def check_wkt_syntax(text):
    """Kontrola vnorenej syntaxe WKT; bez databazy CRS a externych kniznic."""
    tokens = []
    i = 0
    while i < len(text):
        if text[i].isspace():
            i += 1
            continue
        if text[i] == '"':
            i += 1
            while i < len(text):
                if text[i] == '"':
                    if i + 1 < len(text) and text[i + 1] == '"':
                        i += 2
                        continue
                    i += 1
                    break
                i += 1
            else:
                raise ValueError("Neuzavrety retazec vo WKT")
            tokens.append("<string>")
        else:
            match = TOKEN.match(text, i)
            if not match:
                raise ValueError("Nejednoznacny alebo poskodeny WKT")
            tokens.append(match.group())
            i = match.end()

    def value(pos, depth=0):
        if depth > 100 or pos >= len(tokens) or tokens[pos] in "[],()":
            raise ValueError("Neplatna struktura WKT")
        token = tokens[pos]
        pos += 1
        if pos < len(tokens) and tokens[pos] in ("[", "("):
            if not re.fullmatch(r"[A-Za-z_][A-Za-z_0-9]*", token):
                raise ValueError("Neplatny nazov uzla WKT")
            close = "]" if tokens[pos] == "[" else ")"
            pos = value(pos + 1, depth + 1)
            while pos < len(tokens) and tokens[pos] == ",":
                pos = value(pos + 1, depth + 1)
            if pos >= len(tokens) or tokens[pos] != close:
                raise ValueError("Chybajuca ciarka alebo neuzavreta zatvorka vo WKT")
            pos += 1
        return pos

    if len(tokens) < 3 or tokens[0].upper() not in WKT_ROOTS or tokens[1] not in ("[", "("):
        raise ValueError("Zaznam nie je rozpoznany WKT suradnicoveho systemu")
    if value(0) != len(tokens):
        raise ValueError("Nadbytocny obsah za WKT")


def corrected_wkt(data):
    text, separator, padding = data.partition(b"\0")
    if separator and any(padding):
        raise ValueError("Neocakavane data za koncom WKT")
    if not text:
        return None  # Prazdny CRS nie je mozne opravit nahradou uvodzoviek.
    if b"'" not in text:
        return None
    if b'"' in text:
        # Apostrof v standardnom dvojito uzavretom nazve je platny znak.
        check_wkt_syntax(text.decode("utf-8"))
        return None
    if b"''" in text:
        raise ValueError("Nejednoznacne zdvojene apostrofy; automaticka oprava odmietnuta")
    corrected = text.replace(b"'", b'"')
    check_wkt_syntax(corrected.decode("utf-8"))
    return corrected + separator + padding


def inspect(path):
    """Precita iba hlavicky VLR/EVLR, nikdy nedekomprimuje body."""
    patches = []
    records = empty = 0
    with path.open("rb") as file:
        stat = os.fstat(file.fileno())
        header = file.read(375)
        if len(header) < 227 or header[:4] != b"LASF":
            raise ValueError("Neplatna alebo neuplna LAS hlavicka")
        major, minor = header[24:26]
        if major != 1 or minor > 4:
            raise ValueError("Podporovane su LAS verzie 1.0 az 1.4")
        header_size, points, count = struct.unpack_from("<HII", header, 94)
        minimum_header = 375 if minor == 4 else 235 if minor == 3 else 227
        if not minimum_header <= header_size <= points <= stat.st_size:
            raise ValueError("Neplatne posuny v LAS hlavicke")
        areas = [(header_size, count, False, points)]
        if minor == 4:
            start, count = struct.unpack_from("<QI", header, 235)
            if count and not points <= start <= stat.st_size:
                raise ValueError("Neplatny zaciatok EVLR")
            areas.append((start, count, True, stat.st_size))
        for offset, count, extended, limit in areas:
            record_size = 60 if extended else 54
            if count and (offset > limit or count > (limit - offset) // record_size):
                raise ValueError("Neplatny pocet VLR/EVLR")
            for _ in range(count):
                file.seek(offset)
                record = file.read(record_size)
                if len(record) != record_size:
                    raise ValueError("Neuplna hlavicka VLR/EVLR")
                _, user, record_id, length, _ = struct.unpack(
                    "<H16sHQ32s" if extended else "<H16sHH32s", record)
                start = offset + record_size
                end = start + length
                if end > limit:
                    raise ValueError("VLR/EVLR presahuje svoj datovy priestor")
                if record_id == 2112 and user.rstrip(b"\0") in (b"LASF_Projection", b"liblas"):
                    if length > MAX_WKT_BYTES:
                        raise ValueError("WKT je prilis velky na automaticku opravu")
                    data = file.read(length)
                    if len(data) != length:
                        raise ValueError("Neuplny WKT zaznam")
                    records += 1
                    empty += not data.rstrip(b"\0")
                    fixed = corrected_wkt(data)
                    if fixed is not None:
                        patches.append(Patch(start, data, fixed))
                offset = end
    return stat, patches, records, empty


def repair(path, original_stat, patches):
    """Predpokladany hash opravy porovna s celym skutocnym vyslednym suborom."""
    with path.open("r+b") as file:
        stat = os.fstat(file.fileno())
        if (stat.st_size, stat.st_mtime_ns, stat.st_ino) != (
                original_stat.st_size, original_stat.st_mtime_ns, original_stat.st_ino):
            raise ValueError("Subor sa pocas kontroly zmenil")
        expected = hashlib.sha256()
        offset = 0
        while True:
            chunk = file.read(CHUNK_BYTES)
            if not chunk:
                break
            adjusted = bytearray(chunk)
            for patch in patches:
                left = max(patch.offset, offset)
                right = min(patch.offset + len(patch.original), offset + len(chunk))
                if left < right:
                    begin, end = left - patch.offset, right - patch.offset
                    if chunk[left-offset:right-offset] != patch.original[begin:end]:
                        raise ValueError("WKT sa pocas kontroly zmenil")
                    adjusted[left-offset:right-offset] = patch.corrected[begin:end]
            expected.update(adjusted)
            offset += len(chunk)
        if offset != stat.st_size or os.fstat(file.fileno()).st_mtime_ns != stat.st_mtime_ns:
            raise ValueError("Subor sa pocas kontroly zmenil")
        try:
            for patch in patches:
                file.seek(patch.offset)
                file.write(patch.corrected)
            file.flush()
            os.fsync(file.fileno())
            file.seek(0)
            actual = hashlib.sha256()
            while True:
                chunk = file.read(CHUNK_BYTES)
                if not chunk:
                    break
                actual.update(chunk)
            if os.fstat(file.fileno()).st_size != stat.st_size or actual.digest() != expected.digest():
                raise ValueError("Kontrolny sucet vysledku nesuhlasi")
        except BaseException:
            # Povodne metadáta su iba v RAM. Na disk sa ziadna zaloha neuklada.
            for patch in patches:
                file.seek(patch.offset)
                file.write(patch.original)
            file.flush()
            os.fsync(file.fileno())
            raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent,
                        help="Priecinok s LAZ (predvolene priecinok skriptu; bez podpriecinkov)")
    parser.add_argument("--check", action="store_true", help="Iba kontrola, bez zapisu a hashovania bodov")
    args = parser.parse_args(argv)
    if not args.directory.is_dir():
        parser.error("Priecinok neexistuje")
    files = sorted(p for p in args.directory.iterdir()
                   if p.suffix.lower() == ".laz" and p.is_file() and not p.is_symlink())
    changed = unchanged = errors = 0
    for path in files:
        try:
            stat, patches, records, empty = inspect(path)
            if not patches:
                status = "BEZ WKT" if not records else "PRAZDNY WKT" if empty else "BEZ ZMENY"
                print(f"{status}: {path.name}", flush=True)
                unchanged += 1
                continue
            count = sum(a != b for p in patches for a, b in zip(p.original, p.corrected))
            if args.check:
                print(f"NA OPRAVU: {path.name} ({count} uvodzoviek v {len(patches)} zaznamoch)", flush=True)
            else:
                print(f"OPRAVUJEM A OVERUJEM: {path.name} ({stat.st_size / 1024**3:.2f} GiB)", flush=True)
                repair(path, stat, patches)
                print(f"OPRAVENE: {count} uvodzoviek; vsetky ostatne bajty nezmenene", flush=True)
            changed += 1
        except (OSError, ValueError, struct.error) as error:
            errors += 1
            print(f"CHYBA: {path.name}: {error}", file=sys.stderr, flush=True)
    label = "na opravu" if args.check else "opravene"
    print(f"Spolu: {len(files)}, {label}: {changed}, bez zmeny: {unchanged}, chyby: {errors}.", flush=True)
    return 1 if errors else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("\nPrerusene pouzivatelom.", file=sys.stderr)
        sys.exit(130)
