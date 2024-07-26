#ifndef IBDPARSE_H
#define IBDPARSE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>

#define MAX_NAME_LEN 256
#define STEPSIZE 200000
extern int VCF_FRMT_PRESENT; // flag for checking if VCF has the FORMAT field


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


/* Comp_dt - structure for storing comparison data from command-line inputs
    + ids      : all samples and their .hap column indices from genotype input
    + n_ids    : number of samples from genotype input
    + ref_ids  : reference (background) samples 
    + n_refids : number of reference samples
    + uids     : user-specified samples to compare
    + n_uids   : number of user-specified samples
    + upos     : user-specified sites to compare
    + n_upos   : number of user-specified sites
    + uaf      : user-specified allele frequencies
    + n_uaf    : number of user-specified allele frequecies */ 
typedef struct comp_dt {
    Sampl* ids;
    size_t n_ids;
    Sampl* refids;
    size_t n_refids;
    Sampl* uids;
    size_t n_uids;
    unsigned long* upos;
    size_t n_upos;
    Freq* uaf;
    size_t n_uaf;
} Comp_dt;


/* Parses .indv file in IMPUTE input
   Args: Comp_dt* data - pointer to comparison data
         const char* indv_fn - path to .indv file
   Returns: 0 if file parsed successfully, 1 otherwise */
int read_indv(Comp_dt* data, const char* indv_fn);


/* Parses samples to compare from user-provided file
   Args: Comp_dt* data - pointer to comparison data 
         const char* sample_fn - path to sample file
   Returns: 0 if file parsed successfully, 1 if error/no matches found */
int read_sf(Comp_dt* data, const char* sample_fn);


/* Parses samples to compare from command line inputs
   Args: Comp_dt* data - pointer to comparison data 
         const char* sample_str  - sample names as one string (comma-separated)
   Returns: 0 if parsed successfully, 1 if error/no matches found */
int read_scmd(Comp_dt* data, const char* sample_str);


/* Parses reference (background) samples from user-provided file
   Args: Comp_dt* data - pointer to comparison data
         const char* ref_fn - path to reference file
   Returns: 0 if file parsed successfully, 1 if error/no matches found */
int read_rf(Comp_dt* data, const char* ref_fn);


/* Parses sites to compare from user-provided file
   Args: Comp_dt* data - pointer to comparison data 
         const char* pos_fn - path to pos file
         const char* chr - only parse sites on this chromosome (if not NULL)
   Returns: 0 if file parsed successfully, 1 otherwise */
int read_pos(Comp_dt* data, const char* pos_fn, const char* chr);


/* Parses allele frequencies from user-provided file
   Args: Comp_dt* data - pointer to comparison data 
         const char* af_fn - path to allele frequency file
         const char* chr - only parse AF from this chromosome (if not NULL)
   Returns: 0 if file parsed successfully, 1 otherwise */
int read_af(Comp_dt* data, const char* af_fn, const char* chr);


/* Checks if site is in list specified by user
   Args: Comp_dt* data - pointer to comparison data
         const size_t pos - position of site
   Returns: 1 if site found in list; 0 otherwise */
int site_in_upos(Comp_dt* data, const size_t pos);


/* Finds a sample in Comp_dt->ids
   Args: Comp_dt* data - pointer to comparison data
         const char* identifier - sample ID of interest
   Returns: pointer to sample in Comp_dt->ids if found; NULL otherwise */
Sampl* find_sample(Comp_dt* data, const char* identifier);


/* Finds allele frequency at a site in Comp_dt->uaf
   Args: Comp_dt* data - pointer to comparison data
         const size_t pos - position of site of interest
   Returns: pointer to Freq in Comp_dt->uaf if found; NULL otherwise */
Freq* fetch_freq(Comp_dt* data, const size_t pos);


/* Calculates frequency of the alternate allele from input genotypes (IMPUTE)
   Args: const char* hap_str - line in .hap file that corresponds to site
         int n_sample - number of genotype samples
   Returns: frequency of the alternate allele */
double find_f_impute(const char* hap_str, int n_samples);


/* Calculates frequency of the alternate allele from input genotypes (VCF)
   Args: unsigned short* alleles - pointer to allele array at site
         int n_sample - number of genotype samples
   Returns: frequency of the alternate allele */
double find_f_vcf(unsigned short* alleles, int n_samples);


/* Parses sample names from VCF header
   Args: const char* header - pointer to header string
         Comp_dt* data - pointer to comparison data
   Returns: 0 if successfully parsed, 1 otherwise */
int vcf_parse_samples(const char* header, Comp_dt* data);


/* Extracts alleles from genotype fields in VCF
   Args: char* gt - pointer to genotype fields as a single string
         int n_samples - number of genotype samples
         unsigned short* alleles - pointer to alleles array at site
   Returns: 0 if successfully parsed, 1 otherwise */
int vcf_parse_gt(char* gt, int n_samples, unsigned short* alleles);


/* Frees memory allocated by Comp_dt object
   Args: Comp_dt* data - pointer to comparison data
   Returns: 0 on success */
int destroy_CD(Comp_dt* data);

#endif /* IBDPARSE_H */



/* UC Santa Cruz (UCSC) Noncommercial License

ACCEPTANCE
In order to get any license under these terms, you must agree to them as both strict obligations and conditions to all your licenses.

COPYRIGHT LICENSE
The licensor grants you a copyright license for the software to do everything you might do with the software that would otherwise infringe the licensor's copyright in it for any permitted purpose.
However, you may only distribute the software according to Distribution License and make changes or new works based on the software according to Changes and New Works License.

DISTRIBUTION LICENSE
The licensor grants you an additional copyright license to distribute copies of the software. Your license to distribute covers distributing the software with changes and new works permitted by Changes and New Works License.

NOTICES
You must ensure that anyone who gets a copy of any part of the software from you also gets a copy of these terms, as well as the following copyright notice:
This software is Copyright ©2020-2022. The Regents of the University of California (“Regents”). All Rights Reserved.

CHANGES AND NEW WORKS LICENSE
The licensor grants you an additional copyright license to make changes and new works based on the software for any permitted purpose.

PATENT LICENSE
The licensor grants you the right to use the software as described in any patent applications or issued patents resulting from UCSC Case Number 2022-808.

NONCOMMERCIAL PURPOSES
Any noncommercial purpose is a permitted purpose.

COMMERCIAL PURPOSES
Contact Innovation Transfer, UC Santa Cruz, innovation@ucsc.edu , https://officeofresearch.ucsc.edu/iatc/ , for any commercial purpose.

PERSONAL USES
Personal use for research, experiment, and testing for the benefit of public knowledge, personal study, private entertainment, hobby projects, amateur pursuits, or religious observance, without any anticipated commercial application, is use for a permitted purpose.

NONCOMMERCIAL ORGANIZATIONS
Use by any charitable organization, educational institution, public research organization, public safety or health organization, environmental protection organization, or government institution is use for a permitted purpose regardless of the source of funding or obligations resulting from the funding.

FAIR USE
You may have "fair use" rights for the software under the law. These terms do not limit them.

NO OTHER RIGHTS
These terms do not allow you to sublicense or transfer any of your licenses to anyone else, or prevent the licensor from granting licenses to anyone else.  These terms do not imply any other licenses.

PATENT DEFENSE
If you make any written claim that the software infringes or contributes to infringement of any patent, all your licenses for the software granted under these terms end immediately. If your company makes such a claim, all your licenses end immediately for work on behalf of your company.

VIOLATIONS
The first time you are notified in writing that you have violated any of these terms, or done anything with the software not covered by your licenses, your licenses can nonetheless continue if you come into full compliance with these terms, and take practical steps to correct past violations, 
within 32 days of receiving notice.  Otherwise, all your licenses end immediately.

NO LIABILITY
As far as the law allows, the software comes as is, without any warranty or condition, and the licensor will not be liable to you for any damages arising out of these terms or the use or nature of the software, under any kind of legal claim.

DEFINITIONS
The "licensor" is Regents, and the "software" is the software the licensor makes available under these terms.
"You" refers to the individual or entity agreeing to these terms.
"Your company" is any legal entity, sole proprietorship, or other kind of organization that you work for, plus all organizations that have control over, are under the control of, or are under common control with that organization.  
"Control" means ownership of substantially all the assets of an entity, or the power to direct its management and policies by vote, contract, or otherwise.  Control can be direct or indirect.
"Your licenses" are all the licenses granted to you for the software under these terms.
"Use" means anything you do with the software requiring one of your licenses. */
