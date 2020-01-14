#!/usr/bin/env python
# -*- coding: utf-8 -*-
import dask.dataframe as ddf
import pandas as pd
import time
from sklearn.neighbors.kde import KernelDensity
from scipy.optimize import curve_fit
import numpy as np

def infer_distribution_from_contig(contacts, K, K0):
    """
    """
    longest_contig_name = contacts.loc[contacts.P1.idxmax()].N1
    inter_contacts = contacts[
            (contacts.N1 == longest_contig_name)
            & (contacts.N2 == longest_contig_name)]
    inter = np.abs(inter_contacts.P1.values - inter_contacts.P2.values)

    kde = KernelDensity(kernel='gaussian', bandwidth=200).fit(inter.reshape(-1, 1))
    f = lambda x: kde.score_samples(x.reshape(-1, 1))

    # distant
    x1 = np.logspace(np.log10(K0), np.log10(K), 500)
    p = lambda x, a, b: a + b * np.log(x)
    param1, cov = curve_fit(p, x1, f(x1))

    # proximal
    # degree = 30
    # x0 = np.logspace(0, np.log10(K0), 500)
    # param0 = np.polyfit(x0, f(x0), degree)

    # P = (lambda x: np.where(
    #         x < K0,
    #         np.poly1d(param0)(x),
    #         np.where(
    #             x < K,
    #             param1[0] + param1[1] * np.log(x),
    #             param1[0] + param1[1] * np.log(K))
    #         ))
    return param1[0], param1[1]

def main():
    import sys
    if len(sys.argv) != 4:
        print('not enough arguments')
        print('usage: python infer-distribution.py foo.mnd K K0')
        return -1

    mnd_filename = sys.argv[1]
    K = int(sys.argv[2])
    K0 = int(sys.argv[3])

    # print('parsing mnd by dask.dataframe.read_csv', time.time())
    df = ddf.read_csv(
        mnd_filename,
        sep=' ',
        header=None,
        names=['N1', 'P1', 'N2', 'P2'],
        usecols=[1, 2, 5, 6],
        engine='c',
    ).compute()
    # reorder index
    # print('reset indexing', time.time())
    df = df.reset_index(drop=True)

    # print('fitting')
    p = infer_distribution_from_contig(df, K=K, K0=K0)
    print(p[0])
    print(p[1])

main()
