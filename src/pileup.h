#ifndef PILEUP_H
#define PILEUP_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "file-io.h"

#define MAX_FIELD_WIDTH (10240)
#define MAX_COV (128)
#define MAX_ID_LEN (256)
#define STEPSIZE 200000

typedef enum { false, true } bool;

typedef struct pul {
  char chr[ MAX_ID_LEN ];
  unsigned int pos;
  char ref;
  unsigned int cov;
  char bases[MAX_COV];
  size_t best_alt_inx;
  unsigned int base_quals[MAX_COV];
  unsigned int map_quals[MAX_COV];
  int strands[MAX_COV];
} Pul;
typedef struct pul* PulP;

typedef struct pu_chr {
  Pul** puls; /* array of Pul* */
  size_t n_puls;        /* number of pileup lines */
  Pul* pul_arr; /* memories where we'll put actually puls */
} Pu_chr;

char revcom_base( const char base );
size_t base_inx( char b );
int best_base_from_pul( PulP pp,
                        unsigned int mqc,
			unsigned int covc, bool weighted );
int rand_good_base_from_pul( PulP pp,
                             unsigned int mqc,
			     unsigned int covc );

int line2pul( char* line, PulP pp );
int base_inx_from_pul( PulP pp,
                       unsigned int mqc,
		       unsigned int covc  );
unsigned int count_base_from_pul(PulP pul, const char base);
Pul* fetch_Pul( const Pu_chr* puc, const size_t pos );
Pu_chr* init_Pu_chr( const char* fn, const char* chr );
int destroy_Pu_chr(Pu_chr* puc);
#endif /* PILEUP_H */



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