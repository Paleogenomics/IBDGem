#ifndef LOADI2_H
#define LOADI2_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_NAME_LEN 256
#define STEPSIZE 200000

/* Sampl - structure for storing sample information
    + name  : name of sample as appear in .indv file 
    + idx   : index of sample in .indv file */
typedef struct sampl {
    char name[MAX_NAME_LEN+1];
    unsigned int idx;
} Sampl;


/* Freq - structure for storing allele frequencies information
    + pos  : position of variant on chromosome
    + f    : alternate allele frequency of variant */
typedef struct freq {
    unsigned long pos;
    double f;
} Freq;


/* Impute2 - structure for storing genotype information parsed from IMPUTE2-formatted files.
    + haps     : stores phased haplotypes (from .hap file)
    + pos      : stores variant positions (from .legend file)
    + ids      : stores variant rsIDs (from .legend file)
    + ref      : stores reference alleles (from .legend file)
    + alt      : stores alternative alleles (from .legend file)
    + upos     : stores user-specified list of sites to compare
    + n_upos   : number of user-specified sites
    + uaf      : stores user-specified allele frequencies    
    + n_uaf    : number of user-specified allele frequencies
    + samples  : stores sample names and associated .hap column indices (from .indv or user-provided)
    + n_sites  : number of variants
    + n_haps   : number of haplotypes */
typedef struct impute2 {
    unsigned short** haps;
    unsigned long* pos;
    char** ids;
    char** ref;
    char** alt;
    unsigned long* upos;
    size_t n_upos;
    Freq* uaf;
    size_t n_uaf;
    Sampl* samples;
    size_t n_sites;
    size_t n_haps;
} Impute2;


/* Initializes an Impute2 object with genotype information from IMPUTE2-formatted input files
   Args: const char* hap_fn     - path to .hap file
         const char* legend_fn  - path to .legend
         const char* indv_fn    - path to .indv file
   Returns: pointer to an Impute2 object */
Impute2* init_I2(const char* hap_fn, const char* legend_fn, const char* indv_fn);


/* Finds the given sample in i2->samples array
   Args: Impute2* i2             - pointer to Impute genotype info
         const char* identifier  - sample ID of interest
   Returns: pointer to sample in i2->samples if found;
            NULL otherwise */
Sampl* find_sample(Impute2* i2, const char* identifier);


/* Finds allele frequency of a site given its position
   Args: Impute2* i2       - pointer to Impute genotype info
         const size_t pos  - position of site of interest
   Returns: pointer to Freq object in i2->ext_af if found;
            NULL otherwise */
Freq* fetch_freq(Impute2* i2, const size_t pos);


/* Calculates frequency of the alternate allele from input genotypes
   Args: Impute2* i2      - pointer to Impute genotype info
         size_t site_idx  - index position of site
   Returns: frequency of the alternate allele */
double find_f(Impute2* i2, size_t site_idx);


/* Parses samples to compare from a user-provided file
   Args: Impute2* i2            - pointer to Impute genotype info 
         const char* sample_fn  - name of file
   Returns: number of samples read from file;
            -1 if error opening file */
int read_sf(Impute2* i2, const char* sample_fn);


/* Parses samples to compare from command line inputs
   Args: Impute2* i2             - pointer to Impute genotype info 
         const char* sample_str  - sample names as one string (comma-separated)
   Returns: number of samples read from string;
            -1 if error parsing string/no match found */
int read_scmd(Impute2* i2, const char* sample_str);


/* Parses allele frequencies from a user-provided file
   Args: Impute2* i2        - pointer to Impute genotype info 
         const char* af_fn  - name of file
         const char* chr    - only parse AF from this chromosome (if not NULL)
   Returns: number of sites for which allele freqs are read successfully;
            -1 if error opening file */
int read_af(Impute2* i2, const char* af_fn, const char* chr);


/* Parses sites to compare from a user-provided file
   Args: Impute2* i2         - pointer to Impute genotype info 
         const char* pos_fn  - name of file
         const char* chr     - only parse sites on this chromosome (if not NULL)
   Returns: number of sites read successfully;
            -1 if error opening file */
int read_pos(Impute2* i2, const char* pos_fn, const char* chr);


/* Checks if site is in list specified by user
   Args: Impute2* i2       - pointer to Impute genotype info
         const size_t pos  - position of site
   Returns: 1 if site found in list; 0 otherwise */
int site_in_upos(Impute2* i2, const size_t pos);


/* Frees memory allocated by Impute2 object
   Args: Impute2* i2  - pointer to Impute2 object
   Returns: 0 on success */
int destroy_I2(Impute2* i2);

#endif /* LOADI2_H */



/* UC Santa Cruz (UCSC) Noncommercial License

Acceptance
In order to get any license under these terms, you must agree to them as both strict obligations and conditions to all your licenses.

Copyright License
The licensor grants you a copyright license for the software to do everything you might do with the software that would otherwise infringe the licensor's copyright in it for any permitted purpose.
However, you may only distribute the software according to Distribution License and make changes or new works based on the software according to Changes and New Works License.

Distribution License
The licensor grants you an additional copyright license to distribute copies of the software. Your license to distribute covers distributing the software with changes and new works permitted by Changes and New Works License.

Notices
You must ensure that anyone who gets a copy of any part of the software from you also gets a copy of these terms, as well as the following copyright notice:
This software is Copyright ©2020-2022. The Regents of the University of California (“Regents”). All Rights Reserved.

Changes and New Works License
The licensor grants you an additional copyright license to make changes and new works based on the software for any permitted purpose.

Patent License
The licensor grants you the right to use the software as described in any patent applications or issued patents resulting from UCSC Case Number 2022-808.

Noncommercial Purposes
Any noncommercial purpose is a permitted purpose.

Commercial Purposes
Contact Innovation Transfer, UC Santa Cruz, innovation@ucsc.edu , https://officeofresearch.ucsc.edu/iatc/ , for any commercial purpose.

Personal Uses
Personal use for research, experiment, and testing for the benefit of public knowledge, personal study, private entertainment, hobby projects, amateur pursuits, or religious observance, without any anticipated commercial application, is use for a permitted purpose.

Noncommercial Organizations
Use by any charitable organization, educational institution, public research organization, public safety or health organization, environmental protection organization, or government institution is use for a permitted purpose regardless of the source of funding or obligations resulting from the funding.

Fair Use
You may have "fair use" rights for the software under the law. These terms do not limit them.

No Other Rights
These terms do not allow you to sublicense or transfer any of your licenses to anyone else, or prevent the licensor from granting licenses to anyone else.  These terms do not imply any other licenses.

Patent Defense
If you make any written claim that the software infringes or contributes to infringement of any patent, all your licenses for the software granted under these terms end immediately. If your company makes such a claim, all your licenses end immediately for work on behalf of your company.

Violations
The first time you are notified in writing that you have violated any of these terms, or done anything with the software not covered by your licenses, your licenses can nonetheless continue if you come into full compliance with these terms, and take practical steps to correct past violations, 
within 32 days of receiving notice.  Otherwise, all your licenses end immediately.

No Liability
As far as the law allows, the software comes as is, without any warranty or condition, and the licensor will not be liable to you for any damages arising out of these terms or the use or nature of the software, under any kind of legal claim.

Definitions
The "licensor" is Regents, and the "software" is the software the licensor makes available under these terms.
"You" refers to the individual or entity agreeing to these terms.
"Your company" is any legal entity, sole proprietorship, or other kind of organization that you work for, plus all organizations that have control over, are under the control of, or are under common control with that organization.  
"Control" means ownership of substantially all the assets of an entity, or the power to direct its management and policies by vote, contract, or otherwise.  Control can be direct or indirect.
"Your licenses" are all the licenses granted to you for the software under these terms.
"Use" means anything you do with the software requiring one of your licenses. */