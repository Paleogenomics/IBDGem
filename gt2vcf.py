#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import numpy as np


def mergeInfo(gt_fn, info_fn):
    df_gt = pd.read_csv(gt_fn, sep='\t', names=['rsID', 'chrom', 'pos', 'A0', 'A1'],
                        dtype={'rsID': str, 'chrom': int, 'pos': str, 'A0': str, 'A1': str})
    df_info = pd.read_csv(info_fn, sep='\t', comment='#',
                          names=['chrom', 'pos', 'ref', 'alt'], dtype=str)
    df_final = pd.merge(df_gt, df_info[['pos', 'ref', 'alt']], on='pos')
    return df_final


def formatGT(A0, A1, ref, alt):
    if A0 == ref and A1 == ref:
        return '0|0'
    elif A0 == alt and A1 == alt:
        return '1|1'
    elif A0 == ref and A1 == alt:
        return '0|1'
    elif A0 == alt and A1 == ref:
        return '1|0'
    
    
def convertToVCF(df, sampleID, out_fn):
   
    out = open(out_fn, 'w+')
    out.write('##fileformat=VCFv4.1\n')
    out.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    out.write('##contig=<ID=1,assembly=b37,length=249250621>\n')
    out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out.write('##INFO=<ID=VT,Number=.,Type=String,Description="indicates what type of variant the line represents">\n')
    out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n".format(sampleID))
    
    default_QUAL = 100
    default_FILTER = 'PASS'
    default_INFO = 'VT=SNP'
    default_FORMAT = 'GT'
    for idx, row in df.iterrows():
        rsID = row['rsID']
        chrom = row['chrom']
        pos = row['pos']
        ref = row['ref']
        alt = row['alt']
        A0 = row['A0']
        A1 = row['A1']
        gt = formatGT(A0, A1, ref, alt)
        
        lineInfo = tuple()
        lineInfo = (chrom, pos, rsID, ref, alt, default_QUAL, default_FILTER, default_INFO, default_FORMAT,) + (gt,)
        out.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(*lineInfo))

    out.close()
    return None


def main():
    gt_fn = sys.argv[1]
    info_fn = sys.argv[2]
    sampleID = sys.argv[3]
    out_fn = sys.argv[4]
    
    df = mergeInfo(gt_fn, info_fn)
    convertToVCF(df, sampleID, out_fn)

if __name__ == "__main__":
    main()