#### IBDGem
The program in this repository are used to compare the sequence data
in an input IMPUTE2-formated file to the data in an input BAM file. They calculate
the probability of the BAM data given various levels of relatedness to
an individual whose genotypes are described in the IMPUTE file.

### ibdgem.c
```
Find 0, 1, or 2 IBD segments between IMPUTE2 genotype files and BAM file info.
It is assumed that input is on a single chromosome.
This chromosome ID is learned by parsing the genotype file.
-H <HAP file>
-L <LEGEND file>
-I <IDNV file>
-P <MPILEUP file> (converted from BAM)
-D <downsample to this fold-coverage depth>
-S <identifier of sample if there are multiple samples in the INDV file>
-v <if set, make output only for sites that have >=1 variant
    allele in genotype file for this sample>
Format of output table is tab-delimited with columns:
Position, REF_ALLELE, ALT_ALLELE, rsID, AF, DP, VCFA0, VCFA1, BAMnREF, BAMnALT, L(IBD0), L(IBD1), L(IBD2)
```

### summarize-comparison.c.pl
```
Summarizes output of IBDGem
Writes the aggregate probabilities across bins of the input file.
These bins can be plotted to see regional trends.
Output format is: POS_START POS_END AGGR_IBD0 AGGR_IBD1 AGGR_IBD2 NUM_SITES
-n <number of sites per bin; default = 100>
-q <Genotype Quality cutoff; default = 40>
-l <Low coverage cutoff for reference genotype; default = 5>
-h <High coverage cutoff for reference genotype; default = 12>
-L <Low coverage cutoff for comparison data; default = 1>
-H <High coverage cutoff for comparison data; default = 4>
-c <input comparison table file>
```
