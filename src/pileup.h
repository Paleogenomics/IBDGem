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
