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