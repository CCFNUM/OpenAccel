#!/usr/bin/env python3
# File       : plot.py
# Created    : Fri Sep 04 2026 17:03:55 (+0200)
# Author     : Fabian Wermelinger
# Description: General plot utility
# Copyright 2026 Fabian Wermelinger. All Rights Reserved.

import os
import numpy as np
import argparse
import itertools
import matplotlib.pyplot as plt

def parse_args(*, partial=False):
    parser = argparse.ArgumentParser(
        description="General plotting utility for iteration output data")
    # yapf: disable
    parser.add_argument('-d', '--data', nargs='+', type=str, help="ASCII data files", required=True)
    parser.add_argument('-o', '--outfile', type=str, default='plot.png', help="Output filename")
    parser.add_argument('-s', '--show', action='store_true', help="Plot interactively")
    parser.add_argument('--abs', action='store_true', help="Use absolute y-values")
    parser.add_argument('--xcol', type=int, nargs=1, default=0, help="x-column index")
    parser.add_argument('--ycol', type=int, nargs="+", default=[1], help="y-column index(es)")
    parser.add_argument('--xlim', type=float, nargs=2, help="x-limits")
    parser.add_argument('--ylim', type=float, nargs=2, help="y-limits")
    parser.add_argument('--xlabel', type=str, default='x', help="x-label")
    parser.add_argument('--ylabel', type=str, default='y', help="y-label")
    parser.add_argument('--xlog', action='store_true', help="Use logscale for x-axis")
    parser.add_argument('--ylog', action='store_true', help="Use logscale for y-axis")
    parser.add_argument('--linewidth', type=float, default=1.4, help="Linewidth to use")
    parser.add_argument('--markersize', type=float, default=3.0, help="Marker size")
    parser.add_argument('--no-marker', action='store_true', help="Do not plot markers")
    parser.add_argument('--envelope', action='store_true', help="Plot envelope for transient runs")
    parser.add_argument('--grid', action='store_true', help="Draw grid lines")
    # yapf: enable
    if partial:
        return parser.parse_known_args()
    else:
        return parser.parse_args()


def main(args):
    if args.show:
        args.outfile = None
    plot(args.data, args.outfile, args=args)


def readASCII(file, *, comments='#', **kwargs):
    """Read ASCII text file.

    Read various ASCII files and strip comments first.

    Parameters
    ----------
    file : str
        ASCII file path
    comments : str
        Comment character used in file
    kwargs : dict
        Keyword arguments passed to numpy.loadtxt

    Returns
    -------
    NumPy array
    """
    buf = []
    with open(file, 'r') as raw:
        for line in raw:
            l = line.strip()
            if l.startswith(comments):
                continue
            buf.append(l)
    return np.loadtxt(buf, comments=None, **kwargs)


def plot(files, output_file, *, args):
    """Plot ASCII Accel iteration output"""
    fig, ax = plt.subplots()

    marker = itertools.cycle(
        ('o', 's', 'd', '^', '+', 'p', 'v', '.', '<', '>'))

    for f in sorted(files):
        header = None
        with open(f, 'r') as fin:
            for line in fin:
                line = line.strip()
                if line[0] == '#':
                    # last comment contains column labels
                    header = line[1:].strip()
                else:
                    break

        assert header is not None
        data = readASCII(f, comments='#')
        assert(data.shape[1] > 1)
        iter = data[:, args.xcol]
        cols = data[:, args.ycol]

        if args.abs:
            cols = np.absolute(cols)

        # masks
        if args.envelope:
            # mask out start and end iterations of a transient residual
            upper_mask = data[:, 1] == 1  # start inner iteration
            lower_mask = np.append(data[1:, 1] == 1, False)  # end inner iteration

        hsplit = header.split()
        labels = [hsplit[i] for i in args.ycol]
        assert cols.shape[1] == len(labels)
        for i in range(len(labels)):
            lab = labels[i]
            if '[' in lab:  # remove units if present
                lab = lab.split('[')[0].strip()
            kwargs = {'linewidth': args.linewidth}
            if not args.no_marker:
                kwargs['marker'] = next(marker)
                kwargs['ms'] = args.markersize

            if args.envelope:
                h, = ax.plot(iter[lower_mask],
                             cols[lower_mask, i],
                             label=lab,
                             **kwargs)
                ax.plot(iter[upper_mask],
                        cols[upper_mask, i],
                        color=h.get_color(),
                        linestyle=h.get_linestyle(),
                        **kwargs)
            else:
                ax.plot(iter, cols, label=lab, **kwargs)

    if args.xlim is not None:
        ax.set_xlim(args.xlim)
    if args.ylim is not None:
        ax.set_ylim(args.ylim)
    if args.xlog:
        ax.set_xscale(r"log")
    if args.ylog:
        ax.set_yscale(r"log")
    ax.set_xlabel(f"{args.xlabel}")
    ax.set_ylabel(f"{args.ylabel}")
    ax.legend(loc="upper right", fontsize=9)
    plt.grid(args.grid)
    if output_file is None:
        plt.show()
    else:
        kwargs = {}
        if output_file.endswith('.png'):
            kwargs['dpi'] = 600
        fig.savefig(output_file,
                    bbox_inches='tight',
                    pad_inches=1 / 72.27,
                    **kwargs)

if __name__ == "__main__":
    args = parse_args()
    main(args)
