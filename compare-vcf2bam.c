#define EPSILON 0.02
/*
Main function that takes in command-line arguments, processes data
from input files, and prints out results to output file
*/
int main(command-line arguments, input filenames) {
    /* get command-line arguments and filenames, check if they're valid;
    initializes Impute2* i2 with .hap, .legend, and .indv files;
    intializes Pul_chr* pp with mpileup file; */

    /* find column corresponding to the given sample, see below */
    size_t sample_idx = find_sample_idx(i2, sample_identifier from command-line); 

    size_t n = i2->num_sites;
    /* prints header to output file; */
    for (i; i < n; i++) {
        char* str = process_site(i2, pp, i, sample_idx);
        /* prints str to output file; */
    }
    return 0;
}
Two dimensional array of long unsigned int, have Max N as user input, and populate the array before processing the sites


/*
Processes each variant and calculates the information needed for output table
Arguments: pointer to Impute2 object, pointer to Pul_chr object, site index, sample index
Returns: a string containing all the fields needed for the output table
*/
What do we need from the vcf, we need genotype, we can write function to get genotype from vcf or impute2
char* process_site(Impute* i2, Pul_chr* pp, int i, size_t sample_idx) {
    size_t pos = i2->pos[i]; /* pos stored as char* in load_i2, will convert to size_t type */
    char* ref = i2->ref_alleles[i];
    char* alt = i2-> alt_alleles[i];
    char* rsid = i2->ids[i];
    double f = find_f(i2, i);
    unsigned short vcf_A0 = *(i2->haps[i]+sample_idx); // Get allele0 from column corresponding to sample
    unsigned short vcf_A1 = *(i2->haps[i]+(sample_idx+1)); // Get allele1 from the next column
    size_t bam_nREF = count_allele_from_pul(pp, pos, ref[0]);
    size_t bam_nALT = count_allele_from_pul(pp, pos, alt[0]);
    double IBD0 = find_likelihood_Dgivenf(bam_nREF, bam_nALT, f);
    double IBD1 = find_pD_IBD1(bam_nREF, bam_nALT, f);
    double IBD2 = find_pDgivenG(vcf_A0, vcf_A1, bam_nREF, bam_nALT);
    join the fields into one tab-delimited string;
    return string;
}


/* 
Find column in haps matrix that corresponds to the given sample
Arguments: pointer to Impute2 object, sample identifier (ex: HG00096)
Returns: index to haps array that corresponds to the sample
*/
int find_sample_idx(Impute2* i2, char* identifier) {
    
    return sample_idx;
}

unsigned long find_N_choose_k(int N, int k) {

    return 
}
/*
Find population frequency of the alt allele for a given variant
Arguments: pointer to Impute2 object, index of the site of interest
Returns: the population freq of the alt allele
*/
double find_f(Impute2* i2, size_t pos_index) {
    return f;
}

/* Find the number of observed REF or ALT alleles at given position
Arguments: pointer to Pul_chr object, position of variant, allele of interest
Returns: allele count
*/
size_t count_allele_from_pul(Pul_chr* pp, size_t pos, char allele) {
    Pul* p = fetch the right Pul* using the fetch_Pul() function with given pos as key;
    size_t cov = p->cov;
    size_t count = 0;
    for (int i; i < cov; i++) {
        if (p->bases[i] == allele) {
            count++;
        }
    }
    return count;
}

/*
Find likelihood of observed data given genotype
Arguments: genotype ref allele, genotype alt allele, ref allele count, alt allele count
Returns: likelihood of data given genotype
*/
double find_LD_given_G(unsigned short A0, unsigned short A1, size_t nREF, size_t nALT) {
    double L = 1;
    if (nREF == 0 && nALT == 0) {
        return 1;
    }
    unsigned long N_choose_k = find_N_choose_k(nREF, nALT);
    if (A0 == 0 && A1 == 0) {
        L = N_choose_k * ((1-EPSILON)**nREF) * (EPSILON**nALT);
    }
    else if (A0 == 1 && A1 == 1) {
        L = N_choose_k * ((1-EPSILON)**nALT) * (EPSILON**nREF);
    }
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        L = N_choose_k * (0.5**nREF) * (0.5**nALT);
    }
    else {
        fprintf(stderr, "Invalid genotype: A0 = %d, A1 = %d\n", A0, A1);
        exit(1);
    }
    return L;
}

/*
Find likelihood of observed data given population frequency of alt allele
Arguments: ref allele count, alt allele count, population freq of alt allele
Returns: likelihood of data given population freq of alt allele
*/
double find_LD_given_f(size_t nREF, size_t nALT, double f) {
    double L = 1;
    if (nREF == 0 && nALT == 0) {
        return 1;
    }
    L = ((1-f)**2) * find_LD_given_G(0, 0, nREF, nALT) +
        (2 * (1-f) * f) * find_LD_given_G(0, 1, nREF, nALT) +
        (f**2) * find_LD_given_G(1, 1, nREF, nALT);
    return L;
}

/*
Find likelihood of observed data under IBD1 model
Arguments: ref allele count, alt allele count, population freq of alt allele
Returns: likelihood of data under IBD1 model
*/
double find_LD_IBD1(unsigned short A0, unsigned short A1, size_t nREF, size_t nALT, double f) {
    double L;
    if (A0 == 0 && A1 == 0) {
        L = find_LD_given_f(nREF, nALT, f/2);
    }
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        L = (0.5 * find_LD_given_f(nREF, nALT, f/2)
            + 0.5 * find_LD_given_f(nREF, nALT, (f/2)+0.5));
    }
    else if (A0 == 1 && A1 == 1) {
        L = find_LD_given_f(nREF, nALT, (f/2)+0.5);
    }
    return L;
}


