#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 10:37:13 2023

@author: remynguyen
"""
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt

mpl.rcParams['figure.dpi'] = 500
plt.rcParams["font.family"] = "Arial"


def load_args():
    '''Retrieves command line arguments'''
    desc = "plot-chrarm-stats.py: Creates LLR(IBD2/IBD0) and LLR(IBD1/IBD0) scatterplots using aggregated\n \
                                  chromosome-arm LLR values from chrarm-stats.py."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--in', action='store', dest='input', required=True, metavar='INPUT',
                        help="chromosome-arm statistics file")
    parser.add_argument('-o', '--out', action='store', dest='output', required=True, metavar='OUTPUT',
                        help="output file prefix")
    args = parser.parse_args()
    return args


def custom_round(x, base):
    '''Rounds a number to the nearest specified base'''
    return base * round(float(x)/base)


def rename_chrarm(df, arm):
    '''Removes chr prefix and adds p or q suffix to chromosome label'''
    for i in range(len(df)):
        label = df.iloc[i,0][3:] + arm
        df.iloc[i,0] = label
    return df


def read_chrarm_stats(fn):
    '''Reads chrarm-stats.py's LLR output to 4 separate dataframes'''
    df = pd.read_table(fn, delimiter='\t', skiprows=1, names=['CHROM', 'p_20', 'q_20', 'p_10', 'q_10'],
                       dtype={0: str, 1: np.float64, 2: np.float64, 3: np.float64, 4: np.float64},
                       float_precision='high')
    p20 = df[['CHROM', 'p_20']]
    q20 = df[['CHROM', 'q_20']]
    p10 = df[['CHROM', 'p_10']]
    q10 = df[['CHROM', 'q_10']]
    for df in [p20, p10]:
        df = rename_chrarm(df, 'p')
    for df in [q20, q10]:
        df = rename_chrarm(df, 'q')
    return p20, q20, p10, q10


def merge_pq(p_df, q_df):
    '''Merges p- and q-arm dataframes, skipping any NaN'''
    merged_labels = []
    merged_llrs = []
    for i in range(22):
        p_llr = p_df.iloc[i,1]
        if not pd.isna(p_df.iloc[i,1]):
            merged_labels.append(p_df.iloc[i,0])
            merged_llrs.append(p_df.iloc[i,1])
        q_llr = q_df.iloc[i,1]
        if not pd.isna(q_df.iloc[i,1]):
            merged_labels.append(q_df.iloc[i,0])
            merged_llrs.append(q_df.iloc[i,1])
    merged_df = pd.DataFrame({'ARM': merged_labels, 'LLR': merged_llrs})
    return merged_df


def make_plot(df, out):
    '''Creates plot from merged dataframe'''
    fig, ax = plt.subplots(figsize=(6,1.8))
    scatter = ax.scatter(np.arange(len(df)), df['LLR'], color='black', lw=0, s=15)
    ax.spines['top'].set_visible(False)
    ax.set_xticks(np.arange(len(df)))
    ax.set_xticklabels(df['ARM'], fontsize=4.5)
    ax.set_xlabel('chromosome arm', fontsize=7)
    ymin_rounded = custom_round(df['LLR'].min(), 5000)
    ymax_rounded = custom_round(df['LLR'].max(), 5000)
    if ymax_rounded < 0:
        ax.set_yticks(np.arange(ymin_rounded-5000, 5000, 5000))
        ax.set_yticklabels(np.arange(ymin_rounded-5000, 5000, 5000), fontsize=4.5)
    else:
        ax.set_yticks(np.arange(ymin_rounded-5000, ymax_rounded+5000, 5000))
        ax.set_yticklabels(np.arange(ymin_rounded-5000, ymax_rounded+5000, 5000), fontsize=4.5)
    ax.set_ylabel('LLR', fontsize=7)
    ax.axhline(y=0, color='black', linestyle='dotted', alpha=0.5, lw=1)
    ax.tick_params(axis='y', direction='in')
    ymin = df['LLR'].min()
    xmin = np.where(df['LLR']==ymin)[0][0]
    ax.scatter(xmin, ymin, color='red', s=15, lw=0.2, zorder=5, edgecolors='black')
    plt.tight_layout()
    plt.savefig(fname=out+'.png', format='png', dpi=500)
    return None


def main():
    args = load_args()
    fn = args.input
    out = args.output
    p20, q20, p10, q10 = read_chrarm_stats(fn)
    df_20 = merge_pq(p20, q20)
    df_10 = merge_pq(p10, q10)
    make_plot(df_20, out+'.IBD2-IBD0.plot')
    make_plot(df_10, out+'.IBD1-IBD0.plot')


if __name__ == "__main__":
    main()
