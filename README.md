### IBDGem
The program and modules in this repository compare genotype information generated from
deep-sequencing or other means to alignment information from a BAM file. They calculate
the probability of the BAM data given various levels of relatedness to an individual 
whose genotypes are described in the IMPUTE files.

### ibdgem.c
```
Finds 0, 1, or 2 IBD segments between IMPUTE genotype files and BAM file info.

Usage: ./ibdgem -H hap_file -L legend_file -I indv_file -P pileup_file [other options...]
-H (required) <HAP file>
-L (required) <LEGEND file>
-I (required) <INDV file>
-P (required) <PILEUP file>
-N (STR)      <name of individual in Pileup (default: PU_ID)>
-S (FILE)     <file containing subset of samples in panel to
               compare the sequencing data against; one line per sample
-M (INT)      <maximum estimated coverage of Pileup data (default: 20)>
-D (FLOAT)    <down-sample to this fold-coverage depth>
-O (STR)      <path to output directory (default: output to current directory)>
-t (INT)      <number of threads (default: 1)>
-e (FLOAT)    <error rate of sequencing platform used (default: 0.02)>
-v            <if set, make output only for sites that have >=1
               alternate alleles in genotype file for this sample>
               
Format of output table is tab-delimited with columns:
POS, REF, ALT, rsID, AF, DP, VCFA0, VCFA1, BAM_NREF, BAM_NALT, LIBD0, LIBD1, LIBD2
```

### aggregate.c
```
"AGGREGATE: Summarizes the output of IBDGem - partitions the genomic region from the input file
into bins containing a fixed number of SNPs and calculates the aggregated likelihoods for each bin.
These bins can be plotted to see regional trends.
Input data must be sorted by position.

Usage: ./aggregate -c input_table [other options...] >output_file
-c (required) <input comparison table from IBDGem>
-n            <number of sites per bin (default: 100)>
-L            <low coverage cutoff for comparison data (default: 1)>
-H            <high coverage cutoff for comparison data (default: 5)>
Format of output table is tab-delimited with columns:
POS_START POS_END AGGR_LIBD0 AGGR_LIBD1 AGGR_LIBD2 NUM_SITES
```

### hiddengem.c
```
"HIDDENGEM: Find most probable path of IBD states across a genomic region using aggregated likelihoods from IBDGEM.
Usage: ./hiddengem --summary/-s input_summary [other options...] >output_file

--summary/-s (required) <file containing aggregated likelihoods from IBDGem>
--p01                   <penalty for changing between states IBD0 to IBD1 (default: 0.001)>
--p02                   <penalty for changing between states IBD0 to IBD2 (default: 0.000001)>
--p12                   <penalty for changing between states IBD1 to IBD2 (default: 0.001)>
Format of output table is tab-delimited with columns:
Bin, IBD0_Score, IBD1_Score, IBD2_Score, Inferred_State
```

### To compile: 
In the IBDGem directory, type
```
make
```
Separate modules can also be compiled individually by typing
<<<<<<< HEAD
=======
```
### IBDGem
The program and modules in this repository compare genotype information generated from
deep-sequencing or other means to alignment information from a BAM file. They calculate
the probability of the BAM data given various levels of relatedness to an individual 
whose genotypes are described in the IMPUTE files.

### ibdgem.c
```
Finds 0, 1, or 2 IBD segments between IMPUTE genotype files and BAM file info.

Usage: ./ibdgem -H hap_file -L legend_file -I indv_file -P pileup_file [other options...]
-H (required) <HAP file>
-L (required) <LEGEND file>
-I (required) <INDV file>
-P (required) <PILEUP file>
-N (STR)      <name of individual in Pileup (default: PU_ID)>
-S (FILE)     <file containing subset of samples in panel to
               compare the sequencing data against; one line per sample
-M (INT)      <maximum estimated coverage of Pileup data (default: 20)>
-D (FLOAT)    <down-sample to this fold-coverage depth>
-O (STR)      <path to output directory (default: output to current directory)>
-t (INT)      <number of threads (default: 1)>
-e (FLOAT)    <error rate of sequencing platform used (default: 0.02)>
-v            <if set, make output only for sites that have >=1
               alternate alleles in genotype file for this sample>
               
Format of output table is tab-delimited with columns:
POS, REF, ALT, rsID, AF, DP, VCFA0, VCFA1, BAM_NREF, BAM_NALT, LIBD0, LIBD1, LIBD2
```

### aggregate.c
```
"AGGREGATE: Summarizes the output of IBDGem - partitions the genomic region from the input file
into bins containing a fixed number of SNPs and calculates the aggregated likelihoods for each bin.
These bins can be plotted to see regional trends.
Input data must be sorted by position.

Usage: ./aggregate -c input_table [other options...] >output_file
-c (required) <input comparison table from IBDGem>
-n            <number of sites per bin (default: 100)>
-L            <low coverage cutoff for comparison data (default: 1)>
-H            <high coverage cutoff for comparison data (default: 5)>
Format of output table is tab-delimited with columns:
POS_START POS_END AGGR_LIBD0 AGGR_LIBD1 AGGR_LIBD2 NUM_SITES
```

### hiddengem.c
```
"HIDDENGEM: Find most probable path of IBD states across a genomic region using aggregated likelihoods from IBDGEM.
Usage: ./hiddengem --summary/-s input_summary [other options...] >output_file

--summary/-s (required) <file containing aggregated likelihoods from IBDGem>
--p01                   <penalty for changing between states IBD0 to IBD1 (default: 0.001)>
--p02                   <penalty for changing between states IBD0 to IBD2 (default: 0.000001)>
--p12                   <penalty for changing between states IBD1 to IBD2 (default: 0.001)>
Format of output table is tab-delimited with columns:
Bin, IBD0_Score, IBD1_Score, IBD2_Score, Inferred_State
```

### To compile: 
In the IBDGem directory, type
```
make
```
Separate modules can also be compiled individually by typing
```
make [module_name]
```
For example, to compile hiddengem.c, type
```
make hiddengem
```

### To remove object files and executables:
In the IBDGem directory, type
```
make clean
```

### Auxiliary files:
#### gt2vcf.py
```
gt2vcf.py: Converts tab-delimited genotype file (rsID, chrom, pos, A0, A1) to VCF format.
Requires tab-delimited info file with the following fields: chrom, pos, ref, alt for each known SNP.
(See example file in supplementary folder)

Usage: python gt2vcf.py input_genotype input_info sampleID output_filename
```
