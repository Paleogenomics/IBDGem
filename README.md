#### Compare VCF to BAM
This programs in this repository are used to compare the sequence data
in an input VCF file to the data in an input BAM file. They calculate
the probability of the BAM data given various levels of relatedness to
an individual whose genotypes are described in the VCF file.

### compare-vcf2bam-v3.pl
```
compare-vcf2bam.pl Find 0, 1, or 2 IBD segments between
VCF file and bam file info.
It is assumed that input is on a single chromosome.
This chromosome ID is learned by parsing the VCF file.
-V <VCF file>
-B <BAM file>
-D <downsample to this fold-coverage depth>
-Q <Genotype quality minimum; default = 40>
-a <Remove first base observations>
-I <identifier of sample if there are multiple samples in VCF file>
-v <if set, make output only for sites that have >=1 variant
    allele in VCF genotype for this sample>
Format of output table is tab-delimited with columns:
Position, REF_ALLELE, ALT_ALLELE, rsID, AF, GenotypeQ, DP, VCFA0, VCFA1, BAMnREF, BAMnALT, P(IBD0), P(IBD1), P(IBD2)
```

### summarize-comparison.pl
```
summarize-comparison.pl Summarizes output of compare-vcf2bam.pl
Writes the aggregate probabilities across bins of the input file.
These bins can be plotted to see regional trends.
-b <bin size; default = 1000000>
-q <Genotype Quality cutoff; default = 40>
-l <Low coverage cutoff for reference genotype; default = 5>
-h <High coverage cutoff for reference genotype; default = 12>
-L <Low coverage cutoff for comparison data; default = 1>
-H <High coverage cutoff for comparison data; default = 4>
-c <input comparison table file>
```
