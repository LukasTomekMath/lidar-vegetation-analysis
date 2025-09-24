from matplotlib.patches import Polygon

import matplotlib.pyplot as plt
import numpy as np

import argparse
import laspy
import re


def load_laz(file_name: str) -> laspy.LasData:
    with laspy.open(f"..\\..\\data_retiled_papco\\{file_name}") as fh:
        print("Points from Header:", fh.header.point_count)
        las: laspy.LasData = fh.read()
    return las


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest="file_name")

    args = parser.parse_args()
    file_name: str = args.file_name

    regExp = re.compile(r"tile_(?P<upper_x>-\d*)_(?P<upper_y>-\d*)\.laz")
    if (regMatch := regExp.match(file_name)) is None:
        print("Invalid file name!")
        exit(-1)

    upper_x = int(regMatch.group("upper_x"))
    upper_y = int(regMatch.group("upper_y")) + 2000
    area_points = [
        # fmt: off
            (upper_x       , upper_y),
            (upper_x + 2000, upper_y),
            (upper_x + 2000, upper_y - 2000),
            (upper_x       , upper_y - 2000),
        # fmt: on
    ]
    laz_area = Polygon(
        area_points,
        closed=True,
        facecolor=(0, 0, 0, 0),
        edgecolor=(0, 0, 0, 1),
    )

    fig, ax = plt.subplots(1, 1)
    ax.add_patch(laz_area)

    plt.xlim(upper_x - 2000, upper_x + 4000)
    plt.ylim(upper_y - 4000, upper_y + 2000)
    plt.gca().set_aspect("equal")

    las = load_laz(file_name)
    points = las.xyz
    point_count = len(points)
    plt.scatter(points[0 : point_count // 5, 0], points[0 : point_count // 5, 1], s=2)

    for idx, (x, y) in enumerate(area_points, 1):
        plt.scatter(x, y, c="red", s=5)
        plt.text(x, y, str(idx), c="red")

    plt.title(file_name)
    plt.show()


if __name__ == "__main__":
    main()
