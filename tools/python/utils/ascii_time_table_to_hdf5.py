#!/usr/bin/env python3
# File       : ascii_time_table_to_hdf5.py
# Created    : Tue Jan 07 2025 11:16:48 (+0100)
# Author     : Fabian Wermelinger
# Description: Convert ASCII columnar data to h5
# Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.
import numpy as np
import h5py as h5
import argparse


def parse_args(*, partial=False):
    parser = argparse.ArgumentParser(
        description=
        "ASCII to HDF5 converter of columnar data (first column must be monotonically increasing time)"
    )
    # yapf: disable
    parser.add_argument('-f', '--file', type=str, help="ASCII data file", required=True)
    parser.add_argument('-g', '--group', type=str, default="dataset", help="HDF5 group name")
    parser.add_argument('-o', '--outfile', type=str, help="Output filename")
    parser.add_argument('-d', '--delimiter', type=str, default=None, help="Character used to separate values")
    parser.add_argument('-s', '--skiprows', type=int, default= 0, help="Skip the first `skiprows` lines")
    # yapf: enable
    if partial:
        return parser.parse_known_args()
    else:
        return parser.parse_args()


def main(args):
    data = np.loadtxt(args.file,
                      delimiter=args.delimiter,
                      skiprows=args.skiprows)
    fname = args.outfile or args.file + '.h5'
    with h5.File(fname, 'w') as f:
        f.create_dataset(args.group, data=data.transpose())


if __name__ == "__main__":
    args = parse_args()
    main(args)
