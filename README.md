# IBDGem
IBDGem is an identity analysis tool designed to work with low-coverage sequencing data. 
The program compares sequence information from a poor sample (such as a forensic or 
ancient specimen) to genotype information from one or more samples generated independently 
via deep sequencing or microarrays. At each biallelic SNP, IBDGem calculates the probability 
of observing the sequencing data given that they come from an individual who has 0, 1, or 2 
identical-by-descent chromosomes with the person providing the genotypes. In other words, 
the program evaluates the likelihood that the genotypes' source individual could also have 
generated the DNA sample of interest.

## Table of Contents
* [Installation](#installation)
* [Usage](#usage)
* [Main program (ibdgem.c)](#main-program-(ibdgem.c))
* [Summarizing results (aggregate.c)](#summarizing-results-(aggregate.c))
* [Estimating genomic IBD states (hiddengem.c)](#estimating-genomic-IBD-states-(hiddengem.c))
* [Auxiliary files](#auxiliary-files)


## Installation
### To compile
After downloading the latest package from the [Release](https://github.com/Paleogenomics/IBDGem/releases) page, extract the source code file.
Then, navigate to the resulting IBDGem directory and type:
```bash
make
```
Add the IBDGem directory to your ```$PATH``` with:
```bash
export PATH=$(pwd):$PATH
```
Separate modules can also be compiled individually by typing:
```bash
make [module_name]
```
For example, to compile the ```hiddengem``` module, type:
```bash
make hiddengem
```
### To remove object files and executables
In the IBDGem main directory, type:
```bash
make clean
```


## Usage
### Input
The main program (ibdgem.c) calculates likelihoods of IBD states at each SNP. Two types of input are required:
1. Pileup file containing sequence information of the unidentified sample.
2. A VCF *or* 3 files in the IMPUTE reference-panel format (with extensions ```.hap```, ```.legend```, and ```.indv```)
   containing genotype data from a single or multiple test individuals.

The Pileup file can be generated from a BAM file with [samtools](https://www.htslib.org/doc/samtools.html):
```bash
samtools mpileup -o [out.pileup] [in.bam]
```

Using IMPUTE format has the advantage of being able to first filter genotypes by various metrics through [vcftools](https://vcftools.github.io/man_latest.html)
before inputting to IBDGem if desired. The program also runs faster on IMPUTE files.
The 3 IMPUTE files can be generated from a VCF file with:
```bash
vcftools  [ --vcf [in.vcf] | --gzvcf [in.vcf.gz] ] \
          --max-alleles 2 --min-alleles 2 \
          --max-missing 1 \
          --out [out-prefix] \
          --IMPUTE
```

If, instead of a VCF or IMPUTE files, you have a microarray-derived genotype file with 4 columns: 
rsID, chrom, allele1, allele2, you can use the ```gt2vcf.py``` script included in this repository to 
first convert the genotype file to VCF format (see the [Auxiliary files](#auxiliary-files) section).

Then, run the program by customizing the general command:
```bash
ibdgem -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]
```
Or:
```bash
ibdgem -V [vcf-file] -P [pileup-file] [other options...]
```

**Notes**:
1. By default, IBDGem will infer allele frequencies using the genotypes of all individuals in
the VCF/IMPUTE files. This, however, can lead to decreased accuracy in likelihood calculation 
if the number of individuals is small (i.e. fewer than 50). Thus, it is recommended in this case that 
the user provides allele frequencies calculated from a larger reference panel (such as the 1000 Genomes)
in a separate input file via the  ```--allele-freqs``` option.
2. In converting between VCF and IMPUTE, because the ```--IMPUTE``` argument in ```vcftools``` requires
phased data, but IBDGem does not need phase information, one can superficially modify the VCF to
change the genotype notation (```A0/A1``` to ```A0|A1```) with the bash command:
```bash
sed "/^##/! s/\//|/g" unphased.vcf > mockphased.vcf
```
The resulting VCF file can then be converted to IMPUTE normally with ```vcftools --IMPUTE``.

### Output
The output of IBDGem is a tab-delimited file with the following fields: 
**POS, REF, ALT, rsID, AF, DP, GT_A0, GT_A1, SQ_NREF, SQ_NALT, LIBD0, LIBD1, LIBD2**.
Each row corresponds to a single SNP, and the last 3 columns correspond to the likelihoods 
of model IBD0, IBD1, and IBD2, respectively. From this file, the log-likelihood ratio (LLR) 
between any two models at any given site can be calculated. Since sites are considered 
independent, the likelihoods can also be aggregated into non-overlapping bins containing 
a fixed number of sites to increase the discrimatory power between models and measure 
regional signals across the genome (see the [Summarizing results (aggregate.c)](#summarizing-results-(aggregate.c)) section).

## Example run
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

**Note**: For comparisons of real forensic samples, using the ```--variable-sites-only/-v``` option is
recommended to exclude uninformative sites.


## Main program (ibdgem.c)
Below is the full description of the main program's options:
```
IBDGem: Compares low-coverage sequencing data from an unknown sample to known genotypes
        from a reference individual/panel and calculates the likelihood that the samples
        share 0, 1, or 2 IBD chromosomes.

Usage: ./ibdgem -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]
       OR ./ibdgem -V [vcf-file] -P [pileup-file] [other options...]
-V, --vcf  FILE                 VCF file (required if using VCF)
-H, --hap  FILE                 HAP file (required if using IMPUTE)
-L, --legend  FILE              LEGEND file (required if using IMPUTE)
-I, --indv  FILE                INDV file (required if using IMPUTE)
-P, --pileup  FILE              PILEUP file (required)
-N, --pileup-name  STR          Name of individual in Pileup (default: UNKWN)
-A, --allele-freqs  FILE        File containing allele frequencies from an external panel;
                                   must be sorted & whitespace-delimited with columns CHROM, POS, AF;
                                   use in conjunction with -c if includes multiple chromosomes
                                   (default: calculate AF from input genotypes)
-S, --sample-list  FILE         File containing subset of samples to compare the
                                   sequencing data against; one line per sample
                                   (default: compare against all samples in genotype panel)
-s, --sample  STR               Sample(s) to compare the sequencing data against; comma-separated
                                   without spaces if more than one (e.g. sample1,sample2,etc.)
-p, --positions  FILE           List of sites to compare; can be in position list format with 2 columns
                                   CHROM, POS (1-based coordinates) or BED format (0-based coordinates);
                                   use in conjunction with -c if includes multiple chromosomes
                                   (default: perform comparison at all sites)
-q, --min-qual  FLOAT           Genotype quality minimum when using VCF input (default: no minimum)
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


## Summarizing results (aggregate.c)
The ```aggregate``` module summarizes likelihood results generated by the main program.
It takes in the file that *ibdgem* produces and multiplies the likelihoods for each IBD
state over a fixed number of SNPs across the genomic region.
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


## Estimating genomic IBD states (hiddengem.c)
The ```hiddengem``` module estimates the most likely path of IBD states across genomic windows
set by the ```aggregate``` module. It takes in the output of ```aggregate.c``` and returns the
most likely IBD state at each bin/window. The total IBD proportions can then be used to infer
relatedness between samples. 
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


## Auxiliary files:
The ```bin``` directory contains an independent Python script for converting tab-delimited
genotype reports to VCF format, which can then be further converted to IMPUTE if desired.
It requires a reference file with information about the reference and alternate alleles at each
site to fill in the necessary columns in the output VCF.
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