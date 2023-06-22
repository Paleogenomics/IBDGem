### Script for calculating statistics for each chromosome arm from 
### regional aggregated LLRs (e.g. IBDGem's summary output)

import argparse
import numpy as np

## Centromeric regions were determined by choosing the coordinates of the 
## sequence labeled 'centromere' from the 'gap' table on the genome browser, +/- 1Mbp
hg19_centrmrs = {'chr1': [120535435,125535435], 'chr2': [91326172,96326172], 'chr3': [89504855,94504855],
                 'chr4': [48660118,53660118], 'chr5': [45405642,50405642], 'chr6': [57830167,62830167],
                 'chr7': [57054332,62054332], 'chr8': [42838888,47838888], 'chr9': [46367680,51367680],
                 'chr10': [38254936,43254936], 'chr11': [50644206,55644206], 'chr12': [33856695,38856695],
                 'chr13': [1,20000001], 'chr14': [1,20000001], 'chr15': [1,21000001],
                 'chr16': [34335802,39335802], 'chr17': [21263007,26263007], 'chr18': [14460899,19460899],
                 'chr19': [23681783,28681783], 'chr20': [25369570,30369570], 'chr21': [1,15288130],
                 'chr22': [1,17000001]}

## Centromeric regions were determined by choosing the coordinates of the 
## longest model sequence from the 'centromeres' table on the genome browser, +/- 1Mbp
hg38_centrmrs = {'chr1': [121503248,125785433], 'chr2': [91188146,95090558], 'chr3': [90553420,94655575],
                 'chr4': [48712062,52743952], 'chr5': [46309185,50591370], 'chr6': [57553889,60829935],
                 'chr7': [57169654,61828235], 'chr8': [43033745,46877266], 'chr9': [42389636,46518559],
                 'chr10': [38936001,42497441], 'chr11': [50090418,55342400], 'chr12': [33835296,38185253],
                 'chr13': [1,18416385], 'chr14': [1,18538660], 'chr15': [1,20725255],
                 'chr16': [35337667,39265670], 'chr17': [22195019,27566634], 'chr18': [14797856,21561440],
                 'chr19': [23908690,28190875], 'chr20': [25608146,29494540], 'chr21': [1,13280945],
                 'chr22': [1,15419455]}


def load_args():
    '''Retrieves command line arguments'''
    desc = "chrarm-stats.py: Calculates LLR statistics aggregated over whole chromosome arms,\n \
                             excluding centromeric regions."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-s', '--summaries', action='store', dest='summaries', required=True, metavar='SUMMARIES',
                        help="tab-delimited file with columns CHROM and SUMMARY_FILE_PATH (no header) - \
                              use 'chr' prefix in chromosome name; one line per chromosome/file") 
    parser.add_argument('-b', '--build', action='store', dest='build', required=True, choices=['hg38', 'hg19'],
                        metavar='BUILD', help="human genome reference build used (accepted values: hg19, hg38)")
    parser.add_argument('-o', '--out', action='store', dest='out', required=True, metavar='PATH',
                        help="output file prefix")
    args = parser.parse_args()
    return args


def read_sumlist(fp):
    '''Parses input list of summary files'''
    sum_files = dict()
    f = open(fp, 'r')
    line = f.readline()
    while line:
        c, p = line.rstrip().split()
        sum_files[c] = p
        line = f.readline()
    f.close()
    return sum_files


def resolve(value):
    '''Resolves underflow errors by incrementing to the next representable value after 0'''
    if value == 0:
        return np.nextafter(0,1)
    return value


def get_chrarm_stats(build, sum_files, out_prefix):
    '''Sums aggregated LLRs over each chromosome arm'''
    results = dict()
    if build == 'hg38':
        skip_range = hg38_centrmrs
    elif build == 'hg19':
        skip_range = hg19_centrmrs
    out_fn = out_prefix + ".txt"
    out_file = open(out_fn, 'w')
    out_file.write("CHROM\tparm_IBD2/IBD0\tqarm_IBD2/IBD0\tparm_IBD1/IBD0\tqarm_IBD1/IBD0\n")
    for c, p in sum_files.items():
        fn = p.rstrip().split('/')[-1]
        p20 = q20 = p10 = q10 = 0
        sf = open(p, 'r')
        sf.readline()
        line = sf.readline()
        vals = line.rstrip().split("\t")
        start = np.int64(vals[1])
        end = np.int64(vals[2])
        # if 1st bin already overlaps centromere, set values to NA
        if end > skip_range[c][0]:
            p20 = np.nan
            p10 = np.nan

        # sum LLRs on p arm
        while end < skip_range[c][0]:
            p20 += np.log2(resolve(np.float128(vals[5]))/resolve(np.float128(vals[3])))
            p10 += np.log2(resolve(np.float128(vals[4]))/resolve(np.float128(vals[3])))
            line = sf.readline()
            vals = line.rstrip().split()
            start = np.int64(vals[1])
            end = np.int64(vals[2])
        
        # skip centromeric regions
        while start < skip_range[c][1]:    
            line = sf.readline()
            vals = line.rstrip().split()
            start = np.int64(vals[1])
        
        # sum LLRs on q arm
        while line:
            q20 += np.log2(resolve(np.float128(vals[5]))/resolve(np.float128(vals[3])))
            q10 += np.log2(resolve(np.float128(vals[4]))/resolve(np.float128(vals[3])))
            line = sf.readline()
            vals = line.rstrip().split("\t")
        out_file.write("%s\t%.3e\t%.3e\t%.3e\t%.3e\n" % (c, p20, q20, p10, q10))
        sf.close()
    out_file.close()
    return None


def main():
    args = load_args()
    fp = args.summaries
    build = args.build
    out = args.out
    sum_files = read_sumlist(fp)
    get_chrarm_stats(build, sum_files, out)


if __name__ == "__main__":
    main()
