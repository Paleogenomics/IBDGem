#ifndef IBDMATH_H
#define IBDMATH_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>


/* Populates a 2D array of binomial coefficients, or n-choose-k values,
   given a positive integer n for every value k such that n>=k>=0
   Args: unsigned int n  - a positive integer n
   Returns: pointer to 2D array containing nCk values for every integer k */
unsigned long** init_nCk(unsigned int n);


/* Retrieves a specific nCk value given n and k
   Args: unsigned long** nCk  - pointer to nCk 2D array
         unsigned int n       - a positive integer n
         unsigned int k       - a positive integer k
   Returns: the nCk value */
unsigned long retrieve_nCk(unsigned long** nCk, unsigned int n, unsigned int k);


/* Frees memory allocated by nCk object
   Args: unsigned long** nCk  - pointer to nCk array
   Returns: 0 on success */
int destroy_nCk(unsigned long** nCk, unsigned int n);


/* Calculates the probability of observing alignment data
   given genotype data at a single site
   Args: unsigned short A0   - genotyped allele on haplotype 1
         unsigned short A1   - genotyped allele on haplotype 2 
         unsigned int n_ref  - number of aligned reference (0) alleles
         unsigned int n_alt  - number of aligned alternate (1) alleles 
   Returns: probability of seeing the BAM data given the VCF data */ 
double find_pDgG(unsigned long** nCk, double epsilon, unsigned short A0, unsigned short A1, 
                 unsigned int n_ref, unsigned int n_alt);


/* Calculates the probability of observing alignment data
   given population frequency of the alternate allele 
   Args: double f        - population frequency of alternate allele
         double pD_g_00  - probability of data given genotype is homozygous ref
         double pD_g_01  - probability of data given genotype is heterozygous
         double pD_g_11  - probability of data given genotype is homozygous alt
   Returns: probability of seeing BAM data given some 
            background frequency of alternate allele */
double find_pDgf(double f, double pD_g_00, double pD_g_01, double pD_g_11);


/* Calculates the probability of observing alignment data
   given that a site is IBD on 1 chromosome 
   Args: unsigned short A0  - genotyped allele on haplotype 1
         unsigned short A1  - genotyped allele on haplotype 2 
         double f           - population frequency of alternate allele
         double pD_g_00     - probability of data given genotype is homozygous ref
         double pD_g_01     - probability of data given genotype is heterozygous
         double pD_g_11     - probability of data given genotype is homozygous alt
   Returns: probability of seeing BAM data given IBD1 model at site */
double find_pDgIBD1(unsigned short A0, unsigned short A1, double f,
                    double pD_g_00, double pD_g_01, double pD_g_11);

#endif /* IBDMATH_H */
