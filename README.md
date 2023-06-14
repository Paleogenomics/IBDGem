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
* [Example run](#example-run)
* [Main program (ibdgem.c)](#main-program-ibdgemc)
* [Estimating IBD proportions for relatedness detection (hiddengem.c)](#estimating-ibd-proportions-for-relatedness-detection-hiddengemc)
* [Auxiliary files](#auxiliary-files)


## Installation
### To compile
After downloading the latest package from [Releases](https://github.com/Paleogenomics/IBDGem/releases), extract the source code file.
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

If, instead of a VCF or IMPUTE files, you have a general format genotype file with 4 columns: 
rsID, chrom, allele1, allele2, you can use the ```gt2vcf.py``` script included in this repository to 
first convert the genotype file to VCF format (see the [Auxiliary files](#auxiliary-files) section).

#### Linkage disequilibrium (LD) mode
IBDGem version 2.0 provides an option (```--LD```) to take linkage disequilibrium among alleles into account
when calculating the likelihood of the background model (IBD0). To do this, the program uses genotypes
from a reference set of samples, which consists of either all samples from the VCF/IMPUTE files (default)
or a specific subset of those samples specified via ```--background-list```.
IBDGem will then compare the Pileup data against these background individuals and take the average to be
the likelihood of the data under the IBD0 model, over a genomic segment (determined via ```--window-size```).

The program can be run by customizing the general command:
```bash
ibdgem [--LD] -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]
```
Or:
```bash
ibdgem [--LD] -V [vcf-file] -P [pileup-file] [other options...]
```

**Important Notes**:
1. By default, IBDGem will infer allele frequencies using the genotypes of all individuals in
the VCF/IMPUTE files. This, however, can lead to decreased accuracy in likelihood calculation 
if the number of individuals is small (i.e. fewer than 50). Thus, it is recommended in this case that 
the user provides allele frequencies calculated from a larger reference panel (such as the 1000
Genomes) in a separate file via the  ```--allele-freqs``` option. This is especially important
when running the program under the regular, non-LD mode, where likelihoods of the IBD0 model are
calculated on a per-site basis rather than per-haplotype and are thus more dependent on allele
frequencies.

2. When running IBDGem under LD mode, it is important to make sure that the genotype individuals
that are used as background are *unrelated* to each other and to the Pileup individual. If there is
any relatedness, the IBD0 likelihoods will be inflated and LLR(IBD2/IBD0) will be reduced.
In the case that you have a VCF file with several subject individuals that you want to compare the
Pileup data to and a number of reference individuals that you want to use for background model
calculation, you can explicitly specify the list of subject individuals via ```--sample-list``` and
the background individuals via ```--background-list```. On the other hand, if you have only one subject
individual and multiple reference individuals in the VCF, and the Pileup data is suspected to be from
that one subject individual, you can set the Pileup sample name to be the same as the VCF ID of the
subject individual via ```--pileup-name```, and the program will automatically use all other samples in
the VCF as background, without having to explicitly specify them with ```--background-list```.

3. In converting between VCF and IMPUTE, because the ```--IMPUTE``` argument in ```vcftools``` requires
phased data, but IBDGem does not need phase information, one can superficially modify the VCF to
change the genotype notation (```A0/A1``` to ```A0|A1```) with the bash command:
```bash
sed "/^##/! s/\//|/g" unphased.vcf > mockphased.vcf
```
The resulting VCF file can then be converted to IMPUTE normally with ```vcftools --IMPUTE```.

### Output
IBDGem generates 2 tab-delimited files:
(1) Table file (*.tab.txt) with information about each site in the following fields:
  **CHR, rsID, POS, REF, ALT, AF, DP, SQ_NREF, SQ_NALT, GT_A0, GT_A1, LIBD0, LIBD1, LIBD2**.
    Each row corresponds to a single SNP, and the last 3 columns correspond to the likelihoods 
    of model IBD0, IBD1, and IBD2, respectively (these likelihoods are *NOT* calculated under
    LD mode, even with ```--LD``` option specified).
(2) Summary file (*.summary.txt) with information about each genomic segment in the
    following fields: **SEGMENT, START, END, LIBD0, LIBD1, LIBD2, NUM_SITES**. **START** and
    **END** correspond to the physical coordinates of the first and last SNP in the segment,
    respectively. **NUM_SITES** corresponds to the number of SNPs within the segment over
    which the likelihoods for models IBD0, IBD1, and IBD2 are aggregated, which can be set using
    the ```--window-size``` option when running IBDGem (these likelihoods would be calculated
    under LD mode if ```--LD``` option was specified).
    
From these files, the log-likelihood ratio (LLR) between any two models at any given site/segment
can be calculated. For example, in the special case of determining whether the sequence data derives
from the same individual as the genotype data versus the model of it coming from an unrelated individual,
we simply generate LLRs between the IBD2 and IBD0 models.


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

**Note**: In general, using the ```--variable-sites-only/-v``` option is recommended to exclude
uninformative sites.


## Main program (ibdgem.c)
Below is the full description of the main program's options:
```
IBDGem-2.0: Compares low-coverage sequencing data from an unknown sample to known genotypes
            from a reference individual/panel and calculates the likelihood that the samples
            share 0, 1, or 2 IBD chromosomes.

Usage: ./ibdgem [--LD] -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]
       OR ./ibdgem [--LD] -V [vcf-file] -P [pileup-file] [other options...]
--LD                            Linkage disequilibrium mode ON (default: OFF)
-V, --vcf  FILE                 VCF file (required if using VCF)
-H, --hap  FILE                 HAP file (required if using IMPUTE)
-L, --legend  FILE              LEGEND file (required if using IMPUTE)
-I, --indv  FILE                INDV file (required if using IMPUTE)
-P, --pileup  FILE              PILEUP file (required)
-N, --pileup-name  STR          Name of Pileup sample (default: UNKWN)
-A, --allele-freqs  FILE        File containing allele frequencies from a background panel;
                                   must be sorted & whitespace-delimited with columns CHROM, POS, AF;
                                   use in conjunction with -c if includes multiple chromosomes
                                   (default: calculate AF from input genotypes)
-S, --sample-list  FILE         File containing subset of samples to compare the
                                   sequencing data against; one line per sample
                                   (default: compare against all samples in genotype file)
-s, --sample  STR               Sample(s) to compare the sequencing data against; comma-separated
                                   without spaces if more than one (e.g. sample1,sample2,etc.)
-B, --background-list  FILE     File containing subset of samples to be used as the background panel
                                   for calculating IBD0 in LD mode; one line per sample
                                   (default: use all samples in genotype file as background)
-p, --positions  FILE           List of sites to compare; can be in position list format with 2 columns
                                   CHROM, POS (1-based coordinates) or BED format (0-based coordinates);
                                   use in conjunction with -c if includes multiple chromosomes
                                   (default: perform comparison at all sites)
-q, --min-qual  FLOAT           Genotype quality minimum when using VCF input (default: no minimum)
-M, --max-cov  INT              Maximum estimated coverage of Pileup data (default: 20)
-F, --max-af  FLOAT             Maximum alternate allele frequency (default: 1)
-f, --min-af  FLOAT             Minimum alternate allele frequency (default: 0)
-D, --downsample-cov  FLOAT     Down-sample to this fold-coverage depth
-w, --window-size  INT          Number of sites per genomic segment over which likelihood results
                                   are summarized/aggregated (default: 100)
-O, --out-dir  STR              Path to output directory (default: output to current directory)
-c, --chromosome  STR           Chromosome on which the comparison is done; if not specified,
                                   will assume that all inputs are on one single chromosome
-e, --error-rate  FLOAT         Error rate of sequencing platform (default: 0.02)
-v, --variable-sites-only       If set, make output only for sites that are not
                                   homozygous reference in the genotype file for this sample
-h, --help                      Show this help message and exit

Format of likelihood table is tab-delimited with columns:
CHR, rsID, POS, REF, ALT, AF, DP, SQ_NREF, SQ_NALT, GT_A0, GT_A1, LIBD0, LIBD1, LIBD2

Format of summary file is tab-delimited with columns:
SEGMENT, START, END, LIBD0, LIBD1, LIBD2, NUM_SITES
```


## Estimating IBD proportions for relatedness detection (hiddengem.c)
The ```hiddengem``` module estimates the most likely path of IBD states across genomic regions.
It takes in the summary file generated by IBDGem and returns the most likely IBD state at each segment.
The total IBD proportions can then be used to infer relatedness between samples.
```
HIDDENGEM: Finds most probable path of IBD states across genomic segments.

Usage: ./hiddengem -s [summary-file] [other options...] >[out-file]
--summary, -s  FILE      Summary file from IBDGem likelihood calculation (*.summary.txt) (required)
--p01  FLOAT             Penalty for switching between states IBD0 and IBD1 (default: 1e-3)
--p02  FLOAT             Penalty for switching between states IBD0 and IBD2 (default: 1e-6)
--p12  FLOAT             Penalty for switching between states IBD1 and IBD2 (default: 1e-3)
--help                   Show this help message and exit

Format of output table is tab-delimited with columns:
Segment, IBD0_Score, IBD1_Score, IBD2_Score, Inferred_State
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
