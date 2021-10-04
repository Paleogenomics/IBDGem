#include "pileup.h"

char revcom_base( const char base ) {
  switch (base) {
  case 'A' :
    return 'T';
  case 'C' :
    return 'G';
  case 'G' :
    return 'C';
  case 'T' :
    return 'A';
  default :
    return 'N';
  }
}

size_t base_inx( char b ) {
  char B;
  B = toupper(b);
  switch(B) {
  case 'A' :
    return 0;
  case 'C' :
    return 1;
  case 'G' :
    return 2;
  case 'T' :
    return 3;
  default :
    return 4;
  }
}

inline char inx2base( const size_t inx ) {
  switch (inx) {
  case 0 :
    return 'A';
  case 1 :
    return 'C';
  case 2 :
    return 'G';
  case 3 :
    return 'T';
  default :
    return 'N';
  }
}

/* Counts the bases that pass the quality cutoffs
   Returns the single most populous base, unless that is
   an N. If there is a two way tie, returns a random one
   of these two. If there is a three-way tie, returns
   -1.

   Arguments:
   PulP pp -- data structure containing data from one position of pileup
   unsigned int mqc -- minimum mapping quality; that is, don't consider any base
     with a mapping quality score less than this value
   unsigned int covc -- maximum coverage; i.e., return -1 if position coverage
     is greater than this value
   bool weighted -- true if quality scores should be used to weight each base,
     false otherwise

   -1 signifies no good base
*/
int best_base_from_pul( PulP pp, unsigned int mqc,
        unsigned int covc, bool weighted ) {
  size_t base_counts[5];
  size_t cov, i, this_base_index;
  int rand_inx;
  size_t f_counts, s_counts, t_counts, f_inx, s_inx, t_inx;

  cov = pp->cov; // position coverage

  if ( cov > covc ) {
    return -1;
  }

  if ( cov == 0 ) {
    return -1;
  }

  for( i = 0; i <= 4; i++ ) {
    base_counts[i] = 0;
  }

  for ( i = 0; i < pp->cov; i++ ) {
    if ( pp->map_quals[i] >= mqc ) {
        this_base_index = base_inx(pp->bases[i]);

        // if not using weights, add one to the base count for this base. if
        // using weights, add base and map quality scores to base_counts
        if (weighted) {
          base_counts[this_base_index] += pp->base_quals[i] + pp->map_quals[i];
        }
        else {
          base_counts[this_base_index]++;
        }
    }
  }

  /* Find f_inx & s_inx - the indices of the most populous and
     second most populous bases */
  f_counts = 0;
  s_counts = 0;
  t_counts = 0;
  for( i = 0; i <= 4; i++ ) {
    if ( base_counts[i] >= f_counts ) {
      f_inx = i;
      f_counts = base_counts[f_inx];
    }
  }

  /* No or only N good counts? */
  if ( (f_counts == 0) ||
       (f_inx == 4) ) {
    return -1;
  }

  for( i = 0; i <= 4; i++ ) {
    if ( (i != f_inx) &&
         (base_counts[i] >= s_counts) ) {
      s_inx = i;
      s_counts = base_counts[s_inx];
    }
  }
  for( i = 0; i <= 4; i++ ) {
    if ( (i != f_inx) &&
         (i != s_inx) &&
         (base_counts[i] >= t_counts) ) {
      t_inx = i;
      t_counts = base_counts[t_inx];
    }
  }

  /* Three-way tie? */
  if ( f_counts == t_counts ) {
    return -1;
  }

  /* Single winner? */
  if ( f_counts > s_counts ) {
    return f_inx;
  }

  /* Tie? flip a rand() coin */
  if ( (rand() % 2) == 0 ) {
    return f_inx;
  }
  else {
    return s_inx;
  }
}

/* Makes a random ordering of the bases at this position.
   Goes through the bases in this random ordering and returns the
   index of the first one that passes the map and base/strand
   cutoff criteria. If none do, or if the coverage cutoff is
   exceeded at this position, then -1 is returned
*/
int rand_good_base_from_pul( PulP pp,
                             unsigned int mqc,
			     unsigned int covc ) {
  unsigned int rand_ord[MAX_COV];
  int tmp_inx, rand_inx;
  size_t cov;
  size_t i;
  cov = pp->cov;

  if ( cov > covc ) {
    return -1;
  }

  /* First set up the random order array to be ordered 1..n */
  for( i = 0; i < cov; i++ ) {
    rand_ord[i] = i;
  }

  /* Now, randomize it */
  for( i = 0; i < cov; i++ ) {
    /* Pick a random position */
    rand_inx = (rand() % cov);
    /* Switch current position with the random position */
    tmp_inx = rand_ord[i];
    rand_ord[i] = rand_ord[rand_inx];
    rand_ord[rand_inx] = tmp_inx;
  }

  /* Now, go through the random ordering of bases and return the
     index of the first one that passes the cutoff criteria */
  for( i = 0; i < cov; i++ ) {
    if ( pp->map_quals[rand_ord[i]] >= mqc ) {
      return base_inx( pp->bases[rand_ord[i]] );
    }
  }
  return -1;
}

/* line2pul
   Args: char* line : a pileup line of output (with -s flag to samtools view)
         PulP pp    : a pointer to a Pul to be populated with info
                      for this line
   Returns: 0 is everything is copacetic ; 1 if there is a problem
   Ignores (skips past) indels
*/
int line2pul( char* line, PulP pp ) {
  char raw_base_field[MAX_FIELD_WIDTH];
  char raw_base_quals[MAX_FIELD_WIDTH];
  char raw_map_quals[MAX_FIELD_WIDTH];
  int indel_len, i, field_len, cur_base_num;
  char base_code;

  /* First try to just get the first parts of the line to 
     check the coverage. It might be too high and smash the
     stack! */
  if ( sscanf( line, "%s\t%u\t%c\t%u\t",
               pp->chr,
               &pp->pos,
               &pp->ref,
               &pp->cov ) == 4 ) {
    /* Check to make sure there is not more coverage
       than we can handle */
    if ( pp->cov >= MAX_COV ) {
      return 1;
    }
  }
  else { // couldn't even parse the beginning of this line - give up!
    return 2;
  }

  
  if ( sscanf( line, 
               "%s\t%u\t%c\t%u\t%s\t%s\t%s",
               pp->chr,
               &pp->pos,
               &pp->ref,
               &pp->cov,
               raw_base_field,
               raw_base_quals,
               raw_map_quals ) == 7 ) {
    
    if (pp->cov == 0) { // zero coverge site?
      /* Special line with no real data */
      return 0;
    }
      /* Parse the raw_base_field */
    field_len = strlen( raw_base_field );
    i = 0;
    cur_base_num = 0;
    while( i < field_len ) {
      base_code = raw_base_field[i];
      switch( base_code ) {
      case '.' : // Reference base on forward strand
        pp->bases[cur_base_num]   = pp->ref;
        pp->strands[cur_base_num] = 1;
        cur_base_num++;
        i++;
        break;

      case ',' : // Reference base on reverse strand
        pp->bases[cur_base_num]   = pp->ref;
        pp->strands[cur_base_num] = -1;
        cur_base_num++;
        i++;
        break;

      case 'A' :
        pp->bases[cur_base_num]   = 'A';
        pp->strands[cur_base_num] = 1;
        cur_base_num++;
        i++;
        break;

      case 'a' :
        pp->bases[cur_base_num]   = 'A';
        pp->strands[cur_base_num] = -1;
        cur_base_num++;
        i++;
        break;

      case 'C' :
        pp->bases[cur_base_num]   = 'C';
        pp->strands[cur_base_num] = 1;
        cur_base_num++;
        i++;
        break;

      case 'c' :
        pp->bases[cur_base_num]   = 'C';
        pp->strands[cur_base_num] = -1;
        cur_base_num++;
        i++;
        break;

      case 'G' :
        pp->bases[cur_base_num]   = 'G';
        pp->strands[cur_base_num] = 1;
        cur_base_num++;
        i++;
        break;

      case 'g' :
        pp->bases[cur_base_num]   = 'G';
        pp->strands[cur_base_num] = -1;
        cur_base_num++;
        i++;
        break;

      case 'T' :
        pp->bases[cur_base_num]   = 'T';
        pp->strands[cur_base_num] = 1;
        cur_base_num++;
        i++;
        break;

      case 't' :
        pp->bases[cur_base_num]   = 'T';
        pp->strands[cur_base_num] = -1;
        cur_base_num++;
        i++;
        break;

      case 'N' :
        pp->bases[cur_base_num]   = 'N';
        pp->strands[cur_base_num] = 1;
        cur_base_num++;
        i++;
        break;

      case 'n' :
        pp->bases[cur_base_num]   = 'N';
        pp->strands[cur_base_num] = -1;
        cur_base_num++;
        i++;
        break;

      case '-' : // deletion
        i++;
        indel_len = 0;
        base_code = raw_base_field[i];
        while( isdigit(base_code) ) {
          indel_len *= 10;
          indel_len += base_code - 48;
          i++;
          base_code = raw_base_field[i];
        }
        i += indel_len;
        break;

      case '+' : // insertion
        i++;
        indel_len = 0;
        base_code = raw_base_field[i];
        while( isdigit(base_code) ) {
          indel_len *= 10;
          indel_len += base_code - 48;
          i++;
          base_code = raw_base_field[i];
        }
        i += indel_len;
        break;

      case '$' : // end of read segment
        i++;
        break;

      case '^' : // beginning of read segment
        i += 2; // advance past the map-quality character, too
        break;

      case '*' : // deletion marker, treat it like a base
        pp->bases[cur_base_num] = '*';
        pp->strands[cur_base_num] = 0; // no strand for these
        cur_base_num++;
        i++;
        break;

      default :
        fprintf( stderr, "Cannot parse %c in reads field\n",
                 base_code );
        return 1;
      }
    }

    /* Check to see if we got all the bases we were
       expecting */
    if ( cur_base_num != pp->cov ) {
      fprintf( stderr, 
               "Incorrect number of bases read in: %s\n",
               line );
      return 1;
    }
    
    /* Check length of raw_base_quals & raw_map_quals
       fields */
    if ( (strlen( raw_base_quals ) != cur_base_num) &&
         (strlen( raw_map_quals )  != cur_base_num) ) {
      fprintf( stderr, 
               "Incorrect number of base or map quals in: %s\n",
               line );
      return 1;
    }
    
    /* parse raw_base_quals & raw_map_quals field */
    for( i = 0; i < pp->cov; i++ ) {
      pp->base_quals[i] = (raw_base_quals[i] - 33);
      pp->map_quals[i]  = (raw_map_quals[i] - 33);
    }
    return 0;
  }
  
  else {
    return 1;
  }
}

/* Returns a random index of a valid base from this pileup line
   or -1 if the randomly picked base doesn't pass the quality
   cutoffs
 */

int base_inx_from_pul( PulP pp,
                       unsigned int mqc,
		       unsigned int covc  ) {
  int base_num;
  base_num = (rand() % pp->cov);

  if ( (pp->cov <= covc) &&
       (pp->map_quals[base_num] >= mqc) ) {
    return base_num;
  }
  else {
    return -1;
  }
}

Pul* fetch_Pul( const Pu_chr* puc, const size_t pos ) {
  Pul key;
  Pul** found_pul_p;
  key.pos = pos;
  found_pul_p = bsearch( &key,
			 puc->puls,
			 puc->n_puls,
			 sizeof( Pul* ),
			 cmp_Pul );
  if ( found_pul_p == NULL ) {
     return NULL;
  }
  return *found_pul_p;
}

static int cmp_Pul( const void *pul1, const void *pul2 ) {
  Pul* p1 = (Pul*) pul1; // key
  Pul** pp2 = (Pul**) pul2; // pointer to pointer of Pul
  Pul* p2 = *pp2;
  size_t pos1, pos2;
  pos1 = p1->pos;
  pos2 = p2->pos;
  if ( pos1 < pos2 ) {
    return -1;
  }
  if ( pos1 > pos2 ) {
    return 1;
  }
  if ( pos1 == pos2 ) {
    return 0;
  }
  return 2;
}

Pu_chr* init_Pu_chr( const char* fn ) {
  Pu_chr* puc;
  size_t i = 0;
  size_t n;
  int status;
  File_Src* pu_file;
  char pu_str[ MAX_FIELD_WIDTH*2];

  /* How many mpileup lines? Count them! */
  pu_file = init_FS( fn );
  while( get_line_FS( pu_file, pu_str ) != NULL ) {
    i++;
  }
  destroy_FS( pu_file );
  n = i;
  puc = (Pu_chr*)malloc(sizeof(Pu_chr));
  if ( puc == NULL ) {
    return NULL;
  }

  puc->puls = (Pul**)malloc(sizeof(Pul*) * n);
  puc->pul_arr = (Pul*)malloc(sizeof(Pul) * n);
  if ( (puc->puls == NULL) || (puc->pul_arr == NULL) ) {
    return NULL;
  }
  for( i = 0; i < n; i++ ) {
    puc->puls[i] = &puc->pul_arr[i];
  }
  puc->n_puls = n;

  /* Now it's the right size & memories are allocated.
     Let's parse the input file for realsies...*/
  pu_file = init_FS( fn );
  i = 0; // reset index
  while( get_line_FS( pu_file, pu_str ) != NULL ) {
    status = line2pul( pu_str, puc->puls[i] );
    if ( status ) {
      if ( status == 1 ) {
	/* Coverage too high. No big deal */
	;
      }
      if ( status == 2 ) {
	fprintf( stderr, "Problem parsing %s\n", pu_str );
      }
    }
    else {
      i++;
    }
  }
  puc->n_puls = i; // This may be shorter since some lines
  // failed parsing filters, like coverage cutoff
  /* Now, make sure lines are sorted */
  for( i = 0; i < puc->n_puls - 1; i++ ) {
    /* Make sure position of each pul is less than the
       one following, i.e., make sure input was sorted */
    if ( puc->puls[i]->pos > puc->puls[i+1]->pos ) {
      fprintf( stderr, "mpileup lines not sorted!\n" );
      destroy_FS( pu_file );
      return NULL;
    }
  }
  destroy_FS( pu_file );
  return puc;
}

int destroy_Pu_chr(Pu_chr* puc) {
  if (!puc) {
    return 0;
  }
  free(puc->pul_arr);
  free(puc->puls);
  free(puc);
  return 0;
}
