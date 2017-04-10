#!/usr/bin/env python2

import argparse
import os
import yt

def yt_plt(pltPath, annotate, ext):
    """Procedure to save the plots to pltPath.pdf
    """
    ds1= yt.load(pltPath)
    ds1.derived_field_list
    p1 = yt.ProjectionPlot(ds1, "z", "density")
    if annotate:
        p1.annotate_grids()
    p1.save('.'.join((pltPath.rstrip(os.sep), ext)))

def main(args):
    pltPath = args.f
    annotate = args.a
    ext = args.t

    pltPath = os.path.expanduser(pltPath)
    yt_plt(pltPath, annotate, ext)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=main)
    # Args
    parser.add_argument('-f',
                        help='The input directory of a single plot.')
    parser.add_argument('-a',
                        action='store_true',
                        help='Annotate the grids.')
    parser.add_argument('-t',
                        help='The output file extension. e.g. pdf, png, etc.')
    args = parser.parse_args()
    args.func(args)
