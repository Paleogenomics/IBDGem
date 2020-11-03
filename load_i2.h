#ifndef _LOADI2_H_
#define _LOADI2_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define STEPSIZE 200000

/*-----------------------------------------------------------------------------------------------
Impute2 is a structure that stores genotype information parsed from IMPUTE2-formatted files. 
Its members include:
    + haps          : stores phased haplotypes (from .hap file)
    + ids           : stores variant IDs (from .legend file)
    + pos           : stores variant positions (from .legend file)
    + ref_alleles   : stores reference alleles (from .legend file)
    + alt_alleles   : stores alternative alleles (from .legend file)
    + samples       : stores sample names (from .indv file)
    + num_sites     : number of variants
    + num_haps      : number of haplotypes
-----------------------------------------------------------------------------------------------*/
typedef struct impute2 {
    unsigned short** haps;
    char** ids;
    char** pos;
    char** ref_alleles;
    char** alt_alleles;
    char** samples;
    size_t num_sites;
    size_t num_haps;
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