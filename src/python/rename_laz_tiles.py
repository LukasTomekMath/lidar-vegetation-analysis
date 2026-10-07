"""Preview or rename LAZ tiles using header bounds on a 2000 m grid.

Install: python -m pip install "laspy>=2.4,<3"
Preview: python rename_laz_tiles.py "D:\\data"
Rename:  python rename_laz_tiles.py "D:\\data" --apply

Coordinates must be in metres. The grid origin is (0, 0).
Only headers are read; point data is not decompressed. COPC files are skipped.
"""

import argparse
import math
from pathlib import Path
import sys


TILE_SIZE = 2000


def tile_origin(min_x, max_x, min_y, max_y, tolerance=0.01):
    """Find a unique grid cell containing the bounds (within tolerance).

    Ambiguous bounds are rejected instead of guessing from a sparse cloud.
    """
    def axis_origin(lower, upper):
        if not math.isfinite(lower) or not math.isfinite(upper) or lower > upper:
            raise ValueError("Invalid header bounds")
        first = math.ceil((upper - TILE_SIZE - tolerance) / TILE_SIZE)
        last = math.floor((lower + tolerance) / TILE_SIZE)
        if first > last:
            raise ValueError("Bounds cross a 2000 m grid cell boundary")
        if first != last:
            raise ValueError("Bounds do not identify a unique tile")
        return first * TILE_SIZE

    return axis_origin(min_x, max_x), axis_origin(min_y, max_y)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory", type=Path, help="Directory containing LAZ files")
    parser.add_argument("--apply", action="store_true", help="Apply the displayed renames")
    parser.add_argument("--recursive", action="store_true", help="Include subdirectories")
    parser.add_argument("--tolerance", type=float, default=0.01,
                        help="Boundary tolerance in metres (default: 0.01)")
    args = parser.parse_args()
    if not args.directory.is_dir():
        parser.error("The input directory does not exist")
    if not math.isfinite(args.tolerance) or not 0 <= args.tolerance < TILE_SIZE / 2:
        parser.error("Tolerance must be finite and between 0 and 1000 metres (exclusive)")
    try:
        import laspy
    except ImportError:
        parser.exit(1, 'Install the dependency: python -m pip install "laspy>=2.4,<3"\n')

    entries = args.directory.rglob("*") if args.recursive else args.directory.iterdir()
    files = sorted(p for p in entries if p.is_file() and p.suffix.lower() == ".laz"
                   and not p.name.lower().endswith(".copc.laz"))
    if not files:
        print("No LAZ tiles found (COPC files are excluded).")
        return 0

    plan = []
    errors = []
    targets = {}
    for source in files:
        try:
            with laspy.open(source, read_evlrs=False) as reader:
                header = reader.header
                if header.point_count == 0:
                    raise ValueError("Empty point cloud")
                x, y = tile_origin(header.mins[0], header.maxs[0],
                                   header.mins[1], header.maxs[1], args.tolerance)
            target = source.with_name(f"tile_{x}_{y}.laz")
            key = str(target.resolve()).casefold()
            if key in targets:
                raise ValueError(f"Duplicate target name, also requested by {targets[key]}")
            targets[key] = source
            if target != source and target.exists():
                raise ValueError(f"Target already exists: {target}")
            plan.append((source, target))
            status = "UNCHANGED" if source == target else "RENAME"
            print(f"[{status}] {source} -> {target.name}")
        except Exception as exc:
            errors.append(f"{source}: {exc}")

    if errors:
        for error in errors:
            print(f"[ERROR] {error}", file=sys.stderr)
        print("No files renamed. Resolve all errors and run again.", file=sys.stderr)
        return 1

    changes = [(source, target) for source, target in plan if source != target]
    print(f"\nChecked: {len(plan)}; to rename: {len(changes)}")
    if not args.apply:
        print("Preview only. Add --apply to rename files.")
        return 0

    completed = []
    try:
        for source, target in changes:
            # Hard links provide an atomic no-overwrite operation, unlike POSIX rename.
            # Both paths share a directory, so they are on the same filesystem.
            target.hardlink_to(source)
            try:
                source.unlink()
            except OSError:
                target.unlink()
                raise
            completed.append((source, target))
    except OSError as exc:
        print(f"Rename failed: {exc}", file=sys.stderr)
        print(f"Already renamed: {len(completed)} file(s).", file=sys.stderr)
        for source, target in completed:
            print(f"[DONE] {source} -> {target}", file=sys.stderr)
        return 1
    print(f"Renamed {len(completed)} file(s).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
