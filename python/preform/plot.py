#!/usr/bin/env python3
"""Plot a preform profile (x, y) CSV, mirrored about the axis of symmetry.

Usage: python3 plot.py [profile.csv]    (default: example.csv next to this script)
"""

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_xy(path):
    x, y = [], []
    with open(path, newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        next(reader)  # skip header
        for row in reader:
            x.append(float(row[0]))
            y.append(float(row[1]))
    return np.array(x), np.array(y)


def main(argv):
    path = Path(argv[0]) if argv else Path(__file__).with_name("example.csv")
    x, y = read_xy(path)

    plt.plot(x, y, marker=".", linestyle="None", label="profile")
    plt.plot(-x, y, marker=".", linestyle="None", label="mirrored")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Displacement")
    plt.legend()
    plt.xlim([-40, 40])
    plt.ylim([-140, 77])
    plt.gca().set_aspect("equal", adjustable="box")
    plt.show()


if __name__ == "__main__":
    main(sys.argv[1:])
