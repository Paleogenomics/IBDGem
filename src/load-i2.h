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