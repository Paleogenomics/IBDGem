#!/usr/bin/env python3
# -*- coding: utf-8 -*-

### Script for summarizing HiddenGem statistics from multiple chromosomes

import os
import argparse
import pandas as pd
import numpy as np

def load_args():
    '''Retrieves command line arguments'''
    desc = "sum-hiddengem.py: Summarizes HiddenGem outputs and calculates total proportion of the genome\n \
                              shared IBD0, IBD1, and IBD2 between two samples."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--in', action='store', dest='input', required=True, metavar='INPUT',
                        help="tab-delimited file with columns CHROM and HIDDENGEM_FILE_PATH (no header) - \
                              use 'chr' prefix in chromosome name; one line per chromosome/file")
    parser.add_argument('-o', '--out', action='store', dest='output', required=True, metavar='OUTPUT',
                        help="name of output file")
    args = parser.parse_args()
    return args


def read_input(in_file):
    '''Parses list of HiddenGem output files'''
    hg_files = dict()
    f = open(in_file, 'r')
    line = f.readline()
    while line:
        c, p = line.rstrip().split()
        hg_files[c] = p
        line = f.readline()
    f.close()
    return hg_files


def get_hg_stats(hg_files, out):
    '''Calculates per-chromosome & genome-wide IBD proportions'''
    tot_nsegs = tot_ibd0 = tot_ibd1 = tot_ibd2 = 0
    out_file = open(out, 'w')
    out_file.write("CHROM\tN_SEGMENTS\tN_IBD0\tN_IBD1\tN_IBD2\tFRAC_IBD0\tFRAC_IBD1\tFRAC_IBD2\n")
    for c, p in hg_files.items():
        df = pd.read_table(p, header=0, delimiter='\t', comment="#", usecols=['Inferred_State'], dtype={0: int})
        states = df['Inferred_State'].to_numpy()
        n_segs = len(states)
        n_ibd0 =  (states == 0).sum()     
        n_ibd1 = (states == 1).sum()
        n_ibd2 = (states == 2).sum()
        tot_nsegs += n_segs
        tot_ibd0 += n_ibd0
        tot_ibd1 += n_ibd1
        tot_ibd2 += n_ibd2
        out_file.write("{0}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\t{5:.3f}\t{6:.3f}\t{7:.3f}\n".format(c,
                       n_segs, n_ibd0, n_ibd1, n_ibd2, n_ibd0/n_segs, n_ibd1/n_segs, n_ibd2/n_segs))
    pct_ibd0 = (tot_ibd0/tot_nsegs)*100
    pct_ibd1 = (tot_ibd1/tot_nsegs)*100
    pct_ibd2 = (tot_ibd2/tot_nsegs)*100
    out_file.write("# Total segments = {0:d}\n".format(tot_nsegs))
    out_file.write("# Total IBD0 (%) = {0:.3f}\n".format(pct_ibd0))
    out_file.write("# Total IBD1 (%) = {0:.3f}\n".format(pct_ibd1))
    out_file.write("# Total IBD2 (%) = {0:.3f}\n".format(pct_ibd2))
    out_file.close()
    return None


def main():
    args = load_args()
    p = args.input
    out = args.output
    hg_files = read_input(p)
    get_hg_stats(hg_files, out)
    

if __name__ == "__main__":
    main()


