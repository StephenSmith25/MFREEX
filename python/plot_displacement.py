#!/usr/bin/python

import matplotlib.pyplot as plt
import sys,getopt


def main(_dir):

    print _dir

    plt.plot([1,2,3,4], [1,4,9,16], 'ro')
    plt.axis([0, 6, 0, 20])
    plt.savefig("test_rasterization.pdf", dpi=150)


if __name__ == "__main__":
    main(sys.argv[1:])
