#!/usr/bin/env python3
# File       : diff_fields.py
# Created    : Mon Jan 13 2025 13:04:54 (+0100)
# Author     : Fabian Wermelinger
# Description: Diff Exodus II fields
# Copyright 2025 CCFNUM HSLU T&A. All Rights Reserved.
import argparse
import numpy as np
import vtk.numpy_interface.dataset_adapter as dsa
from paraview.simple import *

import argparse


def parse_args(*, partial=False):
    parser = argparse.ArgumentParser(description="Diff Exodus II fields")
    # yapf: disable
    parser.add_argument('-g', '--gold', type=str, help="Exodus data file (gold)", required=True)
    parser.add_argument('-t', '--test', type=str, help="Exodus data file (test)", required=True)
    parser.add_argument('--time_index_gold', type=int, default=-1, help="Time index (gold)")
    parser.add_argument('--time_index_test', type=int, default=-1, help="Time index (test)")
    parser.add_argument('-f', '--fields', nargs='*', type=str, help="Field names to extract")
    parser.add_argument('-l', '--list', action='store_true', help="List available fields")
    # yapf: enable
    if partial:
        return parser.parse_known_args()
    else:
        return parser.parse_args()


def listPointFields(gold, test):
    print("gold:")
    for f in gold.PointVariables.Available:
        print(f'\t{f}')
    print("test:")
    for f in test.PointVariables.Available:
        print(f'\t{f}')

def setTime(step_index):
    tk = GetTimeKeeper()
    src = GetActiveSource()
    if step_index >= 0:
        assert step_index < len(src.TimestepValues)
    else:
        step_index = -1  # last time index
    tk.Time = src.TimestepValues[step_index]
    Show()
    return tk.Time


def getFieldData(field, time_index, reader):
    SetActiveSource(reader)
    t = setTime(time_index)
    raw = servermanager.Fetch(reader)
    wrap = dsa.WrapDataObject(raw)
    assert field in wrap.PointData
    return t, wrap.PointData[field].Arrays[0]


def main(args):
    # read the data
    gold = ExodusIIReader(FileName=args.gold)
    test = ExodusIIReader(FileName=args.test)
    if args.list:
        listPointFields(gold, test)
        return

    if not args.fields:
        print('No fields specified.\nAvailable:')
        listPointFields(gold, test)
        return
    for field in args.fields:
        t_g, d_g = getFieldData(field, args.time_index_gold, gold)
        t_t, d_t = getFieldData(field, args.time_index_test, test)
        d_e = d_t - d_g
        print(f'\nField: {field}')
        print(f'Time gold: {t_g:e}')
        print(f'Time test: {t_t:e}')
        print(f'absolute error min.:', np.min(np.absolute(d_e)), d_e.shape)
        print(f'absolute error max.:', np.max(np.absolute(d_e)), d_e.shape)
        print(f'absolute error norm:', np.linalg.norm(d_e), d_e.shape)
        print(f'{field} norm (gold):', np.linalg.norm(d_g), d_g.shape)
        print(f'{field} norm (test):', np.linalg.norm(d_t), d_t.shape)


if __name__ == "__main__":
    args = parse_args()
    main(args)
