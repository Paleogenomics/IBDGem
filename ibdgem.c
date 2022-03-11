/*** 
 * IBDGEM - Program for genetic identification from low-coverage sequencing data
 * (c) Remy Nguyen & Ed Green / UC Regents
 * Compares genotype information generated from deep-sequencing or other means
 * to alignment information from a BAM file and calculates the probability that
 * the sample in the BAM file shares 0, 1, or 2 IBD chromosomes with the genotype data.
 *  ***/

#include <math.h>
#include <unistd.h>
#include <float.h>
#include <pthread.h>
#include <time.h>

#include "file-io.h"
#include "load_i2.h"
#include "pileup.h"
#include "nchoosek.h"

#define BASES "ACGT"

static double EPSILON = 0.02;
static unsigned int USER_MAX_COV = 20;
static unsigned int NUM_THREADS = 1;

static char* PU_ID = "PU_ID";
static char* OUT = "./";
static int OPT_S = 0; // flag for comparing against only a subset of the VCF samples
static int OPT_V = 0; // flag for skipping sites that are homozygous REF
static int OPT_D = 0; // down-sampling flag
static unsigned long** N_CHOOSE_K;


/* Structure for storing probabilities of observed data & 
other info (ref/alt allele counts, filter flags, etc.) at each site */
typedef struct prob_sum {
    unsigned int failed_filters : 1;
    unsigned int dp;
    double f;
    size_t n_ref;
    size_t n_alt;
    double pD_g_ibd0_f;
    double pD_g_ibd1_00;
    double pD_g_ibd1_01;
    double pD_g_ibd1_11;
    double pD_g_ibd2_00;
    double pD_g_ibd2_01;
    double pD_g_ibd2_11;
} Prob_sum;


/* Structure for storing coverage distribution info */
typedef struct cov_dist {
    unsigned long dist[MAX_COV];
    double mean_cov;
} Cov_dist;


/* Structure for storing do_compare() arguments */
typedef struct comp_args {
    Impute2* i2;
    Prob_sum** ptab;
    Sample_info* samples;
    int n_samples;
    double cull_p;
    Cov_dist* dist;
} Comp_args;


/* Splits the samples to compare (from .indv or a user-provided file) into
subsets to send to individual threads with size depending on the thread index
   Args: Sample_info* all     - all samples against which to compare the sequencing data
         Sample_info* subset  - samples processed by this thread
         int thread_i         - index of thread
         int batch_size       - average number of samples per thread
         int r                - remainder 
         int include_r        - flag for including remainding samples in the subset */
void make_batch(Sample_info* all, Sample_info* subset, int thread_i,
                int batch_size, int r, int include_r) {
    int idx;
    for (int i = 0; i < batch_size; i++) {
        idx = (batch_size * thread_i) + i;
        subset[i] = all[idx];
    }
    if (include_r) {
        for (int i = 0; i < r; i++) {
            idx = (batch_size * (thread_i+1)) + i;
            subset[batch_size+i] = all[idx];    
        }
    }
}


/* Counts how many times an allele occur in a Pileup line
   Args: Pul* pul           - pointer to Pileup line
         const char allele  - allele of interest
   Returns: number of occurrences of the allele */
size_t count_allele_from_pul(Pul* pul, const char allele) {
    size_t count = 0;
    for (int i = 0; i < pul->cov; i++) {
        if (pul->bases[i] == allele) {
            count++;
        }
    }
    return count;
}


/* Calculates ratio for downsampling alignment data to the target 
   depth of coverage
   Args: Pu_chr* puc       - pointer to Pileup alignment info
         double target_dp  - coverage after downsampling
         Cov_dist* idist   - initial coverage distribution
   Returns: the culling ratio */
double find_cull_p(Pu_chr* puc, double target_dp, Cov_dist* idist) {
    double cull_p = 1;
    size_t total_cov = 0;

    // update initial distribution with coverage from Pileup line
    for (int i = 0; i < puc->n_puls; i++) {
        unsigned int cov = puc->puls[i]->cov;
        if (cov <= USER_MAX_COV) {
            total_cov += cov;
            idist->dist[cov]++;
        }
    }
    double mean_cov = (double)total_cov / puc->n_puls;
    idist->mean_cov = mean_cov;
    if (!OPT_D) {
        return cull_p;
    }
    else if (target_dp > mean_cov) {
        fprintf( stderr, "Observed depth is lower than target depth -D. No culling will be done.\n" );
        return cull_p;
    } 
    cull_p = target_dp / mean_cov;
    return cull_p;
}


/* Determines whether a site is a SNP 
   Args: const char* ref  - reference allele at this site 
         const char* alt  - alternate allele at this site
   Returns: 1 if site is single-polymorphic (not indel);
            0 otherwise */
int is_snp(const char* ref, const char* alt) {
    if ((strstr(BASES, ref) && strlen(ref) == 1) && 
        (strstr(BASES, alt) && strlen(alt) == 1)) {
        return 1;
    }
    return 0;
}


/* Finds the given sample in i2->samples array
   Args: Impute2* i2             - pointer to Impute genotype info
         const char* identifier  - sample ID of interest
   Returns: pointer to sample in i2->samples if found;
            NULL otherwise */
Sample_info* find_sample(Impute2* i2, const char* identifier) {
    for (int i = 0; i < (i2->n_haps)/2; i++) {
        if (strcmp(i2->samples[i].name, identifier) == 0) {
            return &i2->samples[i];
        }
    }
    fprintf( stderr, "Sample %s not found in .indv file\n", identifier );
    return NULL;
}


/* Parses the samples to compare from a user-provided file
   Args: Impute2* i2            - pointer to Impute genotype info 
         Sample_info* samples   - samples from file are stored here
         const char* sample_fn  - name of file
   Returns: number of samples read from file;
            -1 if error opening file */
int read_samples(Impute2* i2, Sample_info* samples, const char* sample_fn) {
    File_Src* sf = init_FS(sample_fn);
    if (!sf) {
        return -1;
    }
    char line[MAX_LINE_LEN+1];
    
    int i = 0;
    while (get_line_FS(sf, line)) {
        if (i == i2->n_haps/2) {
            break;
        }
        line[strlen(line)-1] = '\0';
        Sample_info* sp = find_sample(i2, line);
        if (!sp) {
            fprintf( stderr, "Cannot find sample %s in given genotype panel.\n", samples[i].name );
            continue;
        }
        samples[i] = *sp;
        i++;
    }
    if (i == 0) {
        fprintf( stderr, "Error from read_samples(): No matching samples found in %s.\n", sample_fn );
        destroy_FS(sf);
        exit(1);
    }
    destroy_FS(sf);
    return i;
}


/* Calculates frequency of the alternate allele
   Args: Impute2* i2     - pointer to Impute genotype info
         size_t pos_idx  - index position of site
   Returns: frequency of the alternate allele */
double find_f(Impute2* i2, size_t pos_idx) {
    unsigned short* allele_arr = i2->haps[pos_idx];
    size_t n_alt = 0; 
    for (int i = 0; i < i2->n_haps; i++) {
        if (*(allele_arr+i) == 1) {
            n_alt++;
        }
    } 
    return (double)n_alt / i2->n_haps;
}


/* Performs down-sampling of alignment data at a single site
   Args: unsigned int original_count  - initial base coverage at site 
         double cull_p                - ratio for downsampling 
   Returns: down-sampled base coverage */
size_t cull_dp(unsigned int original_count, double cull_p) {
    if (cull_p == 1.0) {
        return original_count;
    }
    size_t new_count = 0;
    for (int i = 0; i < original_count; i++) {
        if ((rand() / (double)RAND_MAX) < cull_p) {
            new_count++;
        }
    }
    return new_count;
}


/* Calculates the probability of observing alignment data
   given genotype data at a single site
   Args: unsigned short A0  - genotyped allele on haplotype 1
         unsigned short A1  - genotyped allele on haplotype 2 
         size_t n_ref       - number of aligned reference (0) alleles
         size_t n_alt       - number of aligned alternate (1) alleles 
   Returns: probability of seeing the BAM data given the VCF data */ 
double find_pDgG(unsigned short A0, unsigned short A1, size_t n_ref, size_t n_alt) {
    double p = 1.0; 
    
    // probability is 1 if there are no observed bases
    if (n_ref == 0 && n_alt == 0) {
        return p;
    }

    unsigned long N_choose_k = retrieve_nCk(N_CHOOSE_K, n_ref+n_alt, n_ref);
    
    if (A0 == 0 && A1 == 0) {
        p = N_choose_k * pow(1-EPSILON, n_ref) * pow(EPSILON, n_alt);
    }
    else if (A0 == 1 && A1 == 1) {
        p = N_choose_k * pow(1-EPSILON, n_alt) * pow(EPSILON, n_ref);
    }

    // If heterozygous site, prob of drawing either allele is 0.5
    // -> binomial function with p=0.5
    // (n choose k) x p**k x (1-p)**(n-k)
    // here n_ref is arbitrarily chosen to be 'successes'
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        p = N_choose_k * pow(0.5, n_ref) * pow(0.5, n_alt);
    }
    else {
        fprintf( stderr, "Invalid genotype: A0 = %d, A1 = %d\n", A0, A1 );
        exit(1);
    }

    // avoids multiplying with 0 during likelihood aggregation
    if (p == 0.0) {
        p = DBL_MIN; 
    }
    return p;
}


/* Calculates the probability of observing alignment data
   given population frequency of the alternate allele 
   Args: double f        - population frequency of alternate allele
         double pD_g_00  - probability of data given genotype is homozygous ref
         double pD_g_01  - probability of data given genotype is heterozygous
         double pD_g_11  - probability of data given genotype is homozygous alt
   Returns: probability of seeing BAM data given some 
            background frequency of alternate allele */
double find_pDgf(double f, double pD_g_00, double pD_g_01, double pD_g_11) {
    double p = 1.0;
    
    // only case where these probs are 1.0 is when nREF+nALT = 0
    if (pD_g_00 == 1 || pD_g_01 == 1 || pD_g_11 == 1) {
        return p;
    }

    // Hardy-Weinberg  P(Data | Genotype)
    p = pow(1-f, 2.0)  *  pD_g_00 +
        2 * (1-f) * f  *  pD_g_01 +
        pow(f, 2.0)    *  pD_g_11;

    if (p == 0.0) {
        p = DBL_MIN;
    }
    return p;
}


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
                    double pD_g_00, double pD_g_01, double pD_g_11) {
    double p = 1.0;

    // For each VCF genotype, use the probability of the *underlying genotypes* 
    // of the BAM data to calculate the probability of the BAM data given that
    // VCF genotype under IBD1 model

    // The probabilities of all 3 underlying BAM genotypes ((0,0), (0,1), (1,1))
    // can be found for any given VCF genotype 
    
    // Homozygous REF: P(G_BAM=(0,1) | G_VCF=(0,0)) * P(D | G_BAM=(0,1)) +
    //                 P(G_BAM=(0,0) | G_VCF=(0,0)) * P(D | G_BAM=(0,0)) + 0
    // -> 0 prob of seeing BAM genotype (1,1) given VCF genotype (0,0)
    if (A0 == 0 && A1 == 0) {
        p = (f * pD_g_01) + ((1-f) * pD_g_00);
    }

    // Heterozygous: P(G_BAM=(0,1) | G_VCF=(0,1)) * P(D | G_BAM=(0,1)) +
    //               P(G_BAM=(0,0) | G_VCF=(0,1)) * P(D | G_BAM=(0,0)) +
    //               P(G_BAM=(1,1) | G_VCF=(0,1)) * P(D | G_BAM=(1,1))
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        p = (0.5 *         pD_g_01) +
            (0.5 * (1-f) * pD_g_00) +
            (0.5 *   f   * pD_g_11);
    }

    // Homozygous ALT: P(G_BAM=(0,1) | G_VCF=(1,1)) * P(D | G_BAM=(0,1)) +
    //                 P(G_BAM=(1,1) | G_VCF=(1,1)) * P(D | G_BAM=(1,1)) + 0
    // -> 0 prob of seeing BAM genotype (0,0) given VCF genotype (1,1)
    else if (A0 == 1 && A1 == 1) {
        p = ((1-f) * pD_g_01) + (f * pD_g_11);
    }

    if (p == 0.0) {
        p = DBL_MIN;
    }
    return p;
}


/* At each site, calculates the probabilities of observing alignment data under
different IBD models given each possible genotype and stores those numbers in ptab 
   Args: Prob_sum** ptab   - pointer to probability info
         Impute2* i2       - pointer to Impute genotype info
         Pu_chr* puc       - pointer to Pileup alignment info
         double target_dp  - coverage after downsampling
         Cov_dist* idist   - initial coverage distribution 
   Returns: culling depth probability */
double do_pD_calc(Prob_sum** ptab, Impute2* i2, Pu_chr* puc, double target_dp,
                  Cov_dist* idist) {
    
    int n = i2->n_sites; 
    double cull_p = find_cull_p(puc, target_dp, idist);

    for (int i = 0; i < n; i++) {
        char* ref = i2->ref[i];
        char* alt = i2->alt[i];
        // not SNP? excluded from analysis
        if (!is_snp(ref, alt)) {
            ptab[i]->failed_filters = 1;
            continue;
        }

        unsigned long pos = i2->pos[i];
        double f = find_f(i2, i);
        Pul* pul = fetch_Pul(puc, pos);
        // not found in Pileup? excluded from analysis
        if (!pul) {
            ptab[i]->failed_filters = 1; 
            continue;
        }
           
        unsigned int cov = pul->cov;
        size_t n_ref = count_allele_from_pul(pul, ref[0]);
        size_t n_alt = count_allele_from_pul(pul, alt[0]);
        // number of observed bases higher than max coverage?
        // excluded from analysis
        if (n_ref+n_alt > USER_MAX_COV) {
            ptab[i]->failed_filters = 1;
            continue;
        }
        
        // store these info for later outputting to file
        n_ref = cull_dp(n_ref, cull_p);
        n_alt = cull_dp(n_alt, cull_p);
        double pD_g_00 = find_pDgG(0, 0, n_ref, n_alt);
        double pD_g_01 = find_pDgG(0, 1, n_ref, n_alt);
        double pD_g_11 = find_pDgG(1, 1, n_ref, n_alt);

        ptab[i]->pD_g_ibd0_f  = find_pDgf(f, pD_g_00, pD_g_01, pD_g_11);
        ptab[i]->pD_g_ibd1_00 = find_pDgIBD1(0, 0, f, pD_g_00, pD_g_01, pD_g_11);
        ptab[i]->pD_g_ibd1_01 = find_pDgIBD1(0, 1, f, pD_g_00, pD_g_01, pD_g_11);
        ptab[i]->pD_g_ibd1_11 = find_pDgIBD1(1, 1, f, pD_g_00, pD_g_01, pD_g_11);
        ptab[i]->pD_g_ibd2_00 = pD_g_00;
        ptab[i]->pD_g_ibd2_01 = pD_g_01;
        ptab[i]->pD_g_ibd2_11 = pD_g_11;
        
        ptab[i]->dp = cov;
        ptab[i]->f = f;
        ptab[i]->n_ref = n_ref;
        ptab[i]->n_alt = n_alt;
    }
    return cull_p;
}


/* Performs comparison between alignment data and the genotype data of given samples;
this function is to be run by each thread */
void* do_compare(void* args) {
    
    Comp_args* cargs = (Comp_args*)args;
    Impute2* i2 = cargs->i2;
    Sample_info* samples = cargs->samples;
    Prob_sum** ptab = cargs->ptab;
    int n_samples = cargs->n_samples;
    double cull_p = cargs->cull_p;
    Cov_dist* idist = cargs->dist;
    
    for (int i = 0; i < n_samples; i++) {
        char out_fn[MAX_FN_LEN+1];
        sprintf( out_fn, "%s/%s.%s.txt", OUT, PU_ID, samples[i].name );
        FILE* fp = fopen(out_fn, "w");
        if (!fp) {
            fprintf( stderr, "Error: Failed to open %s.\n", out_fn );
            exit(1);
        }  
        int idx = samples[i].idx;
        double libd0, libd1, libd2;

        // print the original coverage distribution of Pu sites
        fprintf( fp, "# INITAL COVERAGE DISTRIBUTION:\n" );
        fprintf( fp, "# COVERAGE NUM_SITES\n" );
        for (int cov = 0; cov <= USER_MAX_COV; cov++) {
            fprintf( fp, "# %d %lu\n", cov, idist->dist[cov] );
        }
        fprintf( fp, "# MEAN DEPTH = %lf\n", idist->mean_cov );
        fprintf( fp, "# CULL DEPTH RATIO = %lf\n", cull_p );

        // initialize final coverage distribution
        unsigned long fdist[USER_MAX_COV+1];
        for (int i = 0; i <= USER_MAX_COV; i++) {
            fdist[i] = 0;
        }
        size_t final_total_cov = 0;
        size_t final_nsites = 0;

        fprintf( fp, "# POS\tREF\tALT\trsID\tAF\tDP\tVCFA0\tVCFA1\tBAM_NREF\tBAM_NALT\tL(IBD0)\tL(IBD1)\tL(IBD2)\n" );
        for (int i = 0; i < i2->n_sites; i++) {
        
            unsigned short A0 = *(i2->haps[i]+idx);
            unsigned short A1 = *(i2->haps[i]+(idx+1));
            if (OPT_V && A0 == 0 && A1 == 0) {
                continue; // if -v on, skip homozygous REF sites
            }

            if (ptab[i]->failed_filters) {
                continue; // not SNP, no Pu data or too high cov
            }

            libd0 = ptab[i]->pD_g_ibd0_f;
            if (A0 == 0 && A1 == 0) {
                libd1 = ptab[i]->pD_g_ibd1_00;
                libd2 = ptab[i]->pD_g_ibd2_00;
            }
            else if ( (A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0) ) {
                libd1 = ptab[i]->pD_g_ibd1_01;
                libd2 = ptab[i]->pD_g_ibd2_01;
            }
            else if (A0 == 1 && A1 == 1)  {
                libd1 = ptab[i]->pD_g_ibd1_11;
                libd2 = ptab[i]->pD_g_ibd2_11;
            }
                   
            fprintf( fp, "%zu\t%s\t%s\t%s\t%lf\t%u\t%u\t%u\t%zu\t%zu\t%e\t%e\t%e\n",
                     i2->pos[i], i2->ref[i], i2->alt[i], i2->ids[i], ptab[i]->f, ptab[i]->dp, 
                     A0, A1, ptab[i]->n_ref, ptab[i]->n_alt, libd0, libd1, libd2 );         
                    
            fdist[ptab[i]->n_ref + ptab[i]->n_alt]++;
            final_total_cov += (ptab[i]->n_ref + ptab[i]->n_alt);
		    final_nsites++;
        }
    
        // print the final coverage distribution
        fprintf( fp, "# FINAL COVERAGE DISTRIBUTION:\n" );
        fprintf( fp, "# COVERAGE NUM_SITES\n" );
        for (int cov = 0; cov <= USER_MAX_COV; cov++) {
            fprintf( fp, "# %d %lu\n", cov, fdist[cov] );
        }
        fprintf( fp, "# FINAL MEAN DEPTH = %lf\n", (double)final_total_cov / final_nsites );
        fclose(fp);
    }
    pthread_exit(NULL);
}


/* Prints help message */
void help(int code) {
    fprintf( stderr, "IBDGem: Find 0, 1, or 2 IBD segments between IMPUTE genotype files and BAM file info.\n" );
	fprintf( stderr, "-H (required) <HAP file>\n");
    fprintf( stderr, "-L (required) <LEGEND file>\n");
    fprintf( stderr, "-I (required) <INDV file>\n");
	fprintf( stderr, "-P (required) <PILEUP file>\n" );
    fprintf( stderr, "-N (STR)      <name of individual in Pileup (default: PU_ID)>\n" );
	fprintf( stderr, "-S (FILE)     <file containing subset of samples in panel to\n" );
    fprintf( stderr, "               compare the sequencing data against; one line per sample\n" );
    fprintf( stderr, "               (default: compare against all samples in genotype panel)>\n" );
    fprintf( stderr, "-M (INT)      <maximum estimated coverage of Pileup data (default: 20)>\n" );
    fprintf( stderr, "-D (FLOAT)    <down-sample to this fold-coverage depth>\n" );
    fprintf( stderr, "-O (STR)      <path to output directory (default: output to current directory)>\n" );
    fprintf( stderr, "-t (INT)      <number of threads (default: 1)>\n" );
    fprintf( stderr, "-e (FLOAT)    <error rate of sequencing platform used (default: 0.02)>\n" );
	fprintf( stderr, "-v            <if set, make output only for sites that have >=1\n" );
	fprintf( stderr, "               alternate alleles in genotype file for this sample>\n" );
	fprintf( stderr, "Format of output table is tab-delimited with columns:\n" );
	fprintf( stderr, "POS, REF, ALT, rsID, AF, DP, VCFA0, VCFA1, BAM_NREF, BAM_NALT, LIBD0, LIBD1, LIBD2\n" );
    exit(code);
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {

    clock_t start, end;
    double time_elapsed;
    start = clock();

    int option;
    char* opts = ":H:L:I:P:N:S:M:D:O:t:e:v";

    char hap_fn[MAX_FN_LEN];
    char legend_fn[MAX_FN_LEN];
    char indv_fn[MAX_FN_LEN];
    char pu_fn[MAX_FN_LEN];
    char sample_fn[MAX_FN_LEN];
    double target_dp = 0;

    if (argc == 1) {
        help(0);
        exit(0);
    }
    while ( (option = getopt(argc, argv, opts)) != -1 ) {
        switch (option) {
            case 'H':
                strcpy(hap_fn, optarg);
                break;
            case 'L':
                strcpy(legend_fn, optarg);
                break;
            case 'I':
                strcpy(indv_fn, optarg);
                break;
            case 'P':
                strcpy(pu_fn, optarg);
                break;
            case 'N':
                PU_ID = optarg;
                break;
            case 'S':
                OPT_S = 1;
                strcpy(sample_fn, optarg);
                break;
            case 'M':
                USER_MAX_COV = atoi(optarg);
                break;
            case 't':
                NUM_THREADS = atoi(optarg);
                break;
            case 'e':
                EPSILON = atof(optarg);
                break;
            case 'O':
                OUT = optarg;
                break;
            case 'D':
                target_dp = atof(optarg);
                OPT_D = 1;
                break;
            case 'v':
                OPT_V = 1;
                break;
            case ':':
                fprintf(stderr, "Option -%c missing required argument.\n", optopt);
                exit(0);
            case '?':
                if (isprint(optopt)) {
                    fprintf(stderr, "Invalid option -%c.\n", optopt);
                }
                else {
                    fprintf(stderr, "Invalid option character '\\x%x'.\n", optopt);
                }
                break;
            default:
                help(0);
                exit(0);
        }
    }
    for (int i = optind; i < argc; i++) {
        help(0);
        exit(0);
    }

    Impute2* i2 = init_I2(hap_fn, legend_fn, indv_fn);
    Pu_chr* puc = init_Pu_chr(pu_fn);
    if (!i2 || !puc) {
        help(0);
        exit(1);
    }
    
    int n_samples = i2->n_haps/2;
    Sample_info ids_to_compare[n_samples];
    
    // if user provides file with names of samples to compare
    // against, store them into array
    if (OPT_S) { 
        n_samples = read_samples(i2, ids_to_compare, sample_fn);
        if (n_samples == -1) {
            destroy_I2(i2);
            destroy_Pu_chr(puc);
            exit(0);
        }
    }
    // if file not provided, will do comparison for every sample
    else {
        for (int i = 0; i < n_samples; i++) {
            ids_to_compare[i] = i2->samples[i];
        }
    }

    if ( OPT_D && target_dp <= 0 ) {
        fprintf( stderr, "Invalid target depth of coverage (must be > 0).\n" );
        destroy_I2(i2);
        destroy_Pu_chr(puc);
        exit(0);
    }

    Prob_sum* ptab[i2->n_sites];
    for (int i = 0; i < i2->n_sites; i++) {
        ptab[i] = malloc(sizeof(Prob_sum));
        ptab[i]->failed_filters = 0;
    }
    N_CHOOSE_K = init_nCk(USER_MAX_COV);
    Cov_dist* idist = malloc(sizeof(Cov_dist));
    memset(idist->dist, 0, sizeof(unsigned long) * (USER_MAX_COV+1));
    
    double cull_p = do_pD_calc(ptab, i2, puc, target_dp, idist);

    // too many threads specified? reduce them to match number of samples
    if (NUM_THREADS > n_samples) {
        NUM_THREADS = n_samples;
    }

    Comp_args cargs[NUM_THREADS];
    int r = n_samples % NUM_THREADS;
    // every thread will process at least this many samples
    int batch_size = (n_samples - r) / NUM_THREADS;
    // flag for including the remaining samples in the last thread
    int include_r = 0;
    int thr_nsamples = batch_size;

    int rc;
    pthread_t thr[NUM_THREADS];
    Sample_info* batches[NUM_THREADS];

    for (int ti = 0; ti < NUM_THREADS; ti++) {
        // last thread? if so will also process the remaining samples
        if (ti == NUM_THREADS-1) {
            thr_nsamples = batch_size+r;
            include_r = 1;
        }
        batches[ti] = malloc(thr_nsamples * sizeof(Sample_info));
        make_batch(ids_to_compare, batches[ti], ti, batch_size, r, include_r);

        cargs[ti].i2 = i2;
        cargs[ti].ptab = ptab;
        cargs[ti].samples = batches[ti];
        cargs[ti].n_samples = thr_nsamples;
        cargs[ti].cull_p = cull_p;
        cargs[ti].dist = idist;

        if ((rc = pthread_create(&thr[ti], NULL, do_compare, &cargs[ti]))) {
            fprintf(stderr, "Error: pthread_create, rc: %d\n", rc);
            for (int i = 0; i < i2->n_sites; i++) {
                free(ptab[i]);
            }
            free(batches[ti]);
            free(idist);
            destroy_I2(i2);
            destroy_nCk(N_CHOOSE_K, USER_MAX_COV);
            destroy_Pu_chr(puc);
            return EXIT_FAILURE;
        }
    }

    // block until all threads complete
    for (int ti = 0; ti < NUM_THREADS; ti++) {
        pthread_join(thr[ti], NULL);
        free(batches[ti]);
    }
    
    for (int i = 0; i < i2->n_sites; i++) {
        free(ptab[i]);
    }
    free(idist);
    destroy_I2(i2);
    destroy_nCk(N_CHOOSE_K, USER_MAX_COV);
    destroy_Pu_chr(puc);

    end = clock();
    time_elapsed = ( (double)(end-start)  / CLOCKS_PER_SEC ) / 60;
    fprintf( stderr, "Execution time: %f minutes.\n", time_elapsed );
    return EXIT_SUCCESS;
}
