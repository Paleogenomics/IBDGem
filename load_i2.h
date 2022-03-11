#ifndef _LOADI2_H_
#define _LOADI2_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STEPSIZE 200000
#define MAX_NAME_LEN 256


/*
Sample_info - structure for storing sample information
Members include:
    + name  : name of sample as appear in .indv file 
    + idx   : index of sample in .indv file
*/
typedef struct sample_info {
    char name[MAX_NAME_LEN+1];
    unsigned int idx;
} Sample_info;


/*-----------------------------------------------------------------------------------------------
Impute2 - structure for storing genotype information parsed from IMPUTE2-formatted files. 
Members include:
    + haps     : stores phased haplotypes (from .hap file)
    + pos      : stores variant positions (from .legend file)
    + ids      : stores variant rsIDs (from .legend file)
    + ref      : stores reference alleles (from .legend file)
    + alt      : stores alternative alleles (from .legend file)
    + samples  : stores sample names (from .indv file)
    + n_sites  : number of variants
    + n_haps   : number of haplotypes
-----------------------------------------------------------------------------------------------*/
typedef struct impute2 {
    unsigned short** haps;
    unsigned long* pos;
    char** ids;
    char** ref;
    char** alt;
    Sample_info* samples;
    size_t n_sites;
    size_t n_haps;
} Impute2;


/*
init_I2 - Initialize an Impute2 object with genotype information from IMPUTE2-formatted input files
Arguments   : paths to .hap file, .legend file, and .indv file
Returns     : pointer to an Impute2 object
*/
Impute2* init_I2(const char* hap_fn, const char* legend_fn, const char* indv_fn);

/*
destroy_I2 - Free memory allocated by an Impute2 object
Arguments   : pointer to an Impute2 object
Returns     : 0 on success
*/
int destroy_I2(Impute2* i2);

#endif /*_LOADI2_*/