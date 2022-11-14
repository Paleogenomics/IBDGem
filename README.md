# IBDGem
IBDGem is an identity analysis tool designed to work with low-coverage sequencing data. 
The program compares sequence information from a poor sample (such as a forensic or 
ancient specimen) to genotype information from one or more samples generated independently 
via deep sequencing or microarrays. At each biallelic SNP, IBDGem calculates the probability 
of observing the sequencing data given that they come from an individual who has 0, 1, or 2 
identical-by-descent chromosomes with the person providing the genotypes. In other words, 
the program evaluates the likelihood that the genotypes' source individual could also have 
generated the DNA sample of interest.
 
IBDGem takes in a Pileup file containing sequence information of the unidentified sample, 
as well as genotype data from a test individual or panel of individuals, which consists 
of 3 files in the IMPUTE reference-panel format.
The Pileup file can be generated from a BAM file with ```samtools```:
```bash
samtools mpileup -o [out.pileup] [in.bam]
```

The 3 IMPUTE-formatted files can be generated from a VCF file with ```vcftools```:
```bash
vcftools  [ --vcf [in.vcf] | --gzvcf [in.vcf.gz] ] \
          --max-alleles 2 --min-alleles 2 \
          --out [out-prefix] \
          --IMPUTE
```
The above command will produce 3 files with extensions ```.hap```, ```.legend```, and ```.indv```, 
they are to be used as inputs to IBDGem along with the Pileup file.
If, instead of a VCF, you have a microarray-derived genotype file with 4 columns: rsID, chrom, 
allele1, allele2, you can use the ```gt2vcf.py``` script included in this repository to first 
convert the genotype file to VCF format (see Auxiliary files section below).
**Note**: Because the ```--IMPUTE``` argument in ```vcftools``` assumes phased data, but IBDGem 
does not require phased genotypes, one can superficially modify the VCF to change the genotype
notation from ```A0/A1``` to ```A0|A1``` to work with ```--IMPUTE``` with the bash command:
```bash
sed "/^##/! s/\//|/g" unphased.vcf > mockphased.vcf
```

The output of IBDGem is a tab-delimited file with the following fields: 
POS, REF, ALT, rsID, AF, DP, GT_A0, GT_A1, SQ_NREF, SQ_NALT, LIBD0, LIBD1, LIBD2.
Each row corresponds to a single SNP, and the last 3 columns correspond to the likelihoods 
of model IBD0, IBD1, and IBD2, respectively. From this file, the log-likelihood ratio (LLR) 
between any two models at any given site can be calculated. Since sites are considered 
independent, the likelihoods can also be aggregated into non-overlapping bins containing 
a fixed number of sites to increase the discrimatory power between models and measure 
regional signals across the genome (see the ```aggregate``` module).



### ibdgem.c
```
IBDGem: Compares low-coverage sequencing data from an unknown sample to known genotypes
        from a reference individual/panel and calculates the likelihood that the samples
        share 0, 1, or 2 IBD chromosomes.

Usage: ./ibdgem -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]
-H, --hap  FILE                 HAP file (required)
-L, --legend  FILE              LEGEND file (required)
-I, --indv  FILE                INDV file (required)
-P, --pileup  FILE              PILEUP file (required)
-N, --pileup-name  STR          Name of individual in Pileup (default: UNKWN)
-A, --allele-freqs  FILE        File containing allele frequencies from an external panel;
                                   must be sorted & whitespace-delimited with columns CHROM, POS, AF;
                                   use in conjunction with -c if includes multiple chromosomes
                                   (default: calculate AF from input genotypes)
-S, --sample-list  FILE         File containing subset of samples to compare the
                                   sequencing data against; one line per sample
                                   (default: compare against all samples in genotype panel)
-s, --sample  STR               Sample(s) to compare the sequencing data against, comma-separated
                                   if more than one (e.g. sample1,sample2,etc.)
-p, --positions  FILE           List of sites to compare; can be in position list format with 2 columns
                                   CHROM, POS (1-based coordinates) or BED format (0-based coordinates);
                                   use in conjunction with -c if includes multiple chromosomes
                                   (default: perform comparison at all sites)
-M, --max-cov  INT              Maximum estimated coverage of Pileup data (default: 20)
-F, --max-af  FLOAT             Maximum alternate allele frequency (default: 1)
-f, --min-af  FLOAT             Minimum alternate allele frequency (default: 0)
-D, --downsample-cov  FLOAT     Down-sample to this fold-coverage depth
-O, --out-dir  STR              Path to output directory (default: output to current directory)
-c, --chromosome  STR           Chromosome on which the comparison is done; if not specified,
                                   will assume that all inputs are on one single chromosome
-t, --threads  INT              Number of threads; only recommended when number of samples to compare
                                   exceeds 10,000 (default: 1)
-e, --error-rate  FLOAT         Error rate of sequencing platform (default: 0.02)
-v, --variable-sites-only       If set, make output only for sites that have >=1
                                   alternate alleles in genotype file for this sample
-h, --help                      Show this help message and exit

Format of output table is tab-delimited with columns:
POS, REF, ALT, rsID, AF, DP, GT_A0, GT_A1, SQ_NREF, SQ_NALT, LIBD0, LIBD1, LIBD2
```

### aggregate.c
```
AGGREGATE: Summarizes IBDGem output by partitioning the genomic region into bins containing
           a fixed number of SNPs and calculates the aggregated likelihoods in each bin.
           These bins can be plotted to see regional trends. Input data must be sorted by position.

Usage: ./aggregate -c [comparison-file] [other options...] >[out-file]
-c, --comparison-file  FILE      Comparison table from IBDGem (required)
-n, --num-sites  INT             Number of sites per bin (default: 100)
-H, --max-cov  INT               High coverage cutoff for comparison data (default: 5)
-L, --min-cov  INT               Low coverage cutoff for comparison data (default: 1)
-h, --help                       Show this help message and exit

Format of output table is tab-delimited with columns:
POS_START POS_END AGGR_LIBD0 AGGR_LIBD1 AGGR_LIBD2 NUM_SITES
```

### hiddengem.c
```
HIDDENGEM: Finds most probable path of IBD states across genomic bins from IBDGem-derived
           aggregated likelihoods.

Usage: ./hiddengem -s [summary-file] [other options...] >[out-file]
--summary, -s  FILE      File containing aggregated likelihood results from IBDGem (required)
--p01  FLOAT             Penalty for switching between states IBD0 and IBD1 (default: 1e-3)
--p02  FLOAT             Penalty for switching between states IBD0 and IBD2 (default: 1e-6)
--p12  FLOAT             Penalty for switching between states IBD1 and IBD2 (default: 1e-3)
--help                   Show this help message and exit

Format of output table is tab-delimited with columns:
Bin, IBD0_Score, IBD1_Score, IBD2_Score, Inferred_State
```

## Installation
### To compile
In the IBDGem main directory, type:
```bash
make
```
Separate modules can also be compiled individually by typing:
```bash
make [module_name]
```
For example, to compile the ```hiddengem``` program, type:
```bash
make hiddengem
```
### To remove object files and executables
In the IBDGem main directory, type:
```bash
make clean
```

## Running IBDGem
In the ```supplementary``` directory, you will find the ```ibdgem-test``` folder that contains
a test suite of inputs (as well as example outputs) for test-running the program. The IMPUTE
files store genotype information of 10 SNP sites for 3 samples: ```sample1```, ```sample2```,
and ```sample3```. The ```test1.pileup```, ```test2.pileup```, and ```test3.pileup``` files 
contain the  sequence data for those corresponding samples.

To perform a run on these mock files, for example to compare the sequence data in ```test1.pileup```
to the 3-sample panel whose genotypes are described in the test IMPUTEs, type:
```bash
ibdgem -H test.hap -L test.legend -I test.indv -P test1.pileup -N sample1
```
This will generate 3 output tables for the pairwise comparisons of sample1-vs-sample1/sample2/sample3.
You can find these same tables in the ```output``` folder within ```ibdgem-test```.

For comparisons of real forensic samples, using the ```--variable-sites-only/-v``` option is
recommended to exclude uninformative sites. 

## Auxiliary files:
The ```bin``` directory contains an independent Python script for converting tab-delimited
genotype reports to VCF format, which can then be further converted to IMPUTE2 for use with
IBDGem.
#### gt2vcf.py
```
gt2vcf.py: Converts a tab-delimited genotype file with fields rsID, chrom, allele1, allele2
           (in that order) to VCF format.
           Requires reference info file with fields chrom, pos, ref, alt for each known SNP.
           (See example info file in supplementary folder).

Usage: python gt2vcf.py [genotype-file] [info-file] [sampleID] [out-file]
```
**Note**: ```gt2vcf.py``` requires the Python libraries ```pandas``` and ```numpy```.
These libraries can be installed with ```pip```:
```bash
pip install pandas
pip install numpy
```
To install without root access, add the ```--user``` argument to the commands above.

For any question about the program or suggestions, please reach out to: ```rennguye@ucsc.edu```.