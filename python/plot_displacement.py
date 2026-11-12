#!/usr/bin/env python3
"""Placeholder plotting script (test rasterisation to PDF)."""

import sys

import matplotlib.pyplot as plt


def main(args):
    print(args)

    plt.plot([1, 2, 3, 4], [1, 4, 9, 16], "ro")
    plt.axis([0, 6, 0, 20])
    plt.savefig("test_rasterization.pdf", dpi=150)


if __name__ == "__main__":
    main(sys.argv[1:])
