#ifndef __NCHOOSEK__
#define __NCHOOSEK__

#include <stdio.h>

/* Module for populating a 2D array of binomial coefficients, or n-choose-k values, given a positive integer n
for every value k such that n>=k>=0. */

/*
init_nCk - Initialize a data structure that stores all n-choose-k values given n
Arguments   : a positive integer n
Returns     : pointer to 2D array containing n-choose-k values for every integer k (n>=k>=0)
*/
unsigned long** init_nCk(unsigned int n);

/*
retrieve_nCk - Retrieve a specific n-choose-k value given n and k
Arguments   : pointer to n-choose-k 2D array, a positive integer n,  a positive integer k
Returns     : the appropriate n-choose-k value  
*/
unsigned long retrieve_nCk(unsigned long** nCk, unsigned int n, unsigned int k);

/*
destroy_nCk - Free memory allocated by the n-choose-k data structure
Arguments   : pointer to the n-choose-k array, the integer n used to populate the structure
Returns     : 0 on success
*/
int destroy_nCk(unsigned long** nCk, unsigned int n);

#endif /*_NCHOOSEK_*/