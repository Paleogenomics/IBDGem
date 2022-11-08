/*** 
 * IBDGEM - Program for genetic identification from low-coverage sequencing data
 * Compares sequence alignment information from an unknown sample to known genotypes
 * generated through deep-sequencing or other means and evaluates the likelihood that
 * the sample in the sequencing data shares 0, 1, or 2 IBD chromosomes with the
 * individual(s) providing the genotype data.
 *  ***/

#include <unistd.h>
#include <pthread.h>
#include <getopt.h>
#include <time.h>

#include "file-io.h"
#include "load-i2.h"
#include "pileup.h"
#include "ibd-math.h"

#define BASES "ACGT"

static double EPSILON = 0.02;
static unsigned int USER_MAX_COV = 20;
static double USER_MAX_AF = 1;
static double USER_MIN_AF = 0;
static unsigned int NUM_THREADS = 1;
static char* SQ_ID = "UNKWN";

static int OPT_S1 = 0; // flag for comparing a subset of samples from user-specified file
static int OPT_S2 = 0; // flag for comparing a subset of samples supplied through command line
static int OPT_A = 0; // flag for using user-specified allele frequencies
static int OPT_P = 0; // flag for using user-specified sites
static int OPT_V = 0; // flag for skipping sites that are homozygous REF
static int OPT_D = 0; // down-sampling flag

/** Long options table **/
static struct option longopts[] = {
    { "hap",                  required_argument, 0, 'H' },
    { "legend",               required_argument, 0, 'L' },
    { "indv",                 required_argument, 0, 'I' },
    { "pileup",               required_argument, 0, 'P' },
    { "pileup-name",          required_argument, 0, 'N' }, 
    { "allele-freqs",         required_argument, 0, 'A' },
    { "sample-list",          required_argument, 0, 'S' },
    { "sample",               required_argument, 0, 's' },
    { "max-cov",              required_argument, 0, 'M' },
    { "downsample-cov",       required_argument, 0, 'D' },
    { "out-dir",              required_argument, 0, 'O' },
    { "max-af",               required_argument, 0, 'F' },
    { "min-af",               required_argument, 0, 'f' },
    { "positions",            required_argument, 0, 'p' },
    { "chromosome",           required_argument, 0, 'c' },
    { "threads",              required_argument, 0, 't' },
    { "error-rate",           required_argument, 0, 'e' },
    { "variable-sites-only",  no_argument, 0, 'v' },
    { "help",                 no_argument, 0, 'h' },
    { 0, 0, 0, 0}
};

/* Structure for storing probabilities of observed data & 
other info (ref/alt allele counts, filter flags, etc.) at each site */
typedef struct prob_sum {
    unsigned int failed_filters : 1;
    unsigned int dp;
    double f;
    unsigned int n_ref;
    unsigned int n_alt;
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
    char* str_cmd;
    Impute2* i2;
    Prob_sum** ptab;
    Sampl* samples;
    int n_samples;
    double cull_p;
    Cov_dist* dist;
} Comp_args;


/* Splits the samples to compare (from .indv or a user-provided file) into
subsets to send to individual threads with size depending on the thread index
   Args: Sampl* all      - all samples against which to compare the sequencing data
         Sampl* subset   - samples processed by this thread
         int thread_i    - index of thread
         int batch_size  - average number of samples per thread
         int r           - remainder 
         int include_r   - flag for including remainding samples in the subset */
void make_batch(Sampl* all, Sampl* subset, int thread_i,
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


/* Performs down-sampling of alignment data at a single site
   Args: unsigned int original_count  - initial base coverage at site 
         double cull_p                - ratio for downsampling 
   Returns: down-sampled base coverage */
unsigned int cull_dp(unsigned int original_count, double cull_p) {
    if (cull_p == 1.0) {
        return original_count;
    }
    unsigned int new_count = 0;
    for (int i = 0; i < original_count; i++) {
        if ((rand() / (double)RAND_MAX) < cull_p) {
            new_count++;
        }
    }
    return new_count;
}


/* At each site, calculates the probabilities of observing alignment data under
different IBD models given each possible genotype and stores those numbers in ptab 
   Args: unsigned long** nCk  - pointer to the nchoosek matrix
         Prob_sum** ptab      - pointer to probability info
         Impute2* i2          - pointer to Impute genotype info
         Pu_chr* puc          - pointer to Pileup alignment info
         double target_dp     - coverage after downsampling
         Cov_dist* idist      - initial coverage distribution 
   Returns: culling depth probability */
double do_pD_calc(unsigned long** nCk, Prob_sum** ptab, Impute2* i2,
                  Pu_chr* puc, double target_dp, Cov_dist* idist) {
    
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
        if (OPT_P) {
            // not a user-specified site? excluded from analysis
            if (!site_in_upos(i2, pos)) {
                ptab[i]->failed_filters = 1;
                continue;
            }
        }
        Pul* pul = fetch_Pul(puc, pos);
        // not found in Pileup? excluded from analysis
        if (!pul) {
            ptab[i]->failed_filters = 1;
            continue;
        }

        double f;
        if (OPT_A) {
            Freq* fp = fetch_freq(i2, pos);
            if (!fp) {
                // if not found, calculate from genotypes
                f = find_f(i2, i);
            }
            else {
                f = fp->f;
            }
        }
        else {
            f = find_f(i2, i);
        }
        // allele frequency higher or lower than specified?
        // excluded from analysis
        if (f > USER_MAX_AF || f < USER_MIN_AF) {
            ptab[i]->failed_filters = 1;
            continue;
        }
           
        unsigned int cov = pul->cov;
        unsigned int n_ref = count_base_from_pul(pul, ref[0]);
        unsigned int n_alt = count_base_from_pul(pul, alt[0]);
        // number of observed bases higher than max coverage?
        // u know the drill
        if (n_ref+n_alt > USER_MAX_COV) {
            ptab[i]->failed_filters = 1;
            continue;
        }

        // store these info for later outputting to file
        n_ref = cull_dp(n_ref, cull_p);
        n_alt = cull_dp(n_alt, cull_p);
        double pD_g_00 = find_pDgG(nCk, EPSILON, 0, 0, n_ref, n_alt);
        double pD_g_01 = find_pDgG(nCk, EPSILON, 0, 1, n_ref, n_alt);
        double pD_g_11 = find_pDgG(nCk, EPSILON, 1, 1, n_ref, n_alt);

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
    Sampl* samples = cargs->samples;
    Prob_sum** ptab = cargs->ptab;
    int n_samples = cargs->n_samples;
    double cull_p = cargs->cull_p;
    Cov_dist* idist = cargs->dist;
    
    for (int i = 0; i < n_samples; i++) {

        // go to the specified output dir
        char out_fn[MAX_FN_LEN];
        sprintf( out_fn, "%s.%s.txt", SQ_ID, samples[i].name );
        FILE* fp = fopen(out_fn, "w");
        if (!fp) {
            fprintf( stderr, "[::] ERROR in do_compare(): Failed to open %s.\n", out_fn );
            exit(1);
        }

        int idx = samples[i].idx;
        double libd0, libd1, libd2;

        // print command line arguments
        fprintf( fp, "# Entered command: %s\n\n", cargs->str_cmd );
        
        // print the original coverage distribution of Pu sites
        fprintf( fp, "# INITAL COVERAGE DISTRIBUTION:\n" );
        fprintf( fp, "# COVERAGE N_SITES\n" );
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

        fprintf( fp, "# POS\tREF\tALT\trsID\tAF\tDP\tGT_A0\tGT_A1\tSQ_NREF\tSQ_NALT\tLIBD0\tLIBD1\tLIBD2\n" );
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
                   
            fprintf( fp, "%zu\t%s\t%s\t%s\t%lf\t%u\t%u\t%u\t%u\t%u\t%e\t%e\t%e\n",
                     i2->pos[i], i2->ref[i], i2->alt[i], i2->ids[i], ptab[i]->f, ptab[i]->dp, 
                     A0, A1, ptab[i]->n_ref, ptab[i]->n_alt, libd0, libd1, libd2 );         
                    
            fdist[ptab[i]->n_ref + ptab[i]->n_alt]++;
            final_total_cov += (ptab[i]->n_ref + ptab[i]->n_alt);
		    final_nsites++;
        }
    
        // print the final coverage distribution
        fprintf( fp, "# FINAL COVERAGE DISTRIBUTION:\n" );
        fprintf( fp, "# COVERAGE N_SITES\n" );
        for (int cov = 0; cov <= USER_MAX_COV; cov++) {
            fprintf( fp, "# %d %lu\n", cov, fdist[cov] );
        }
        fprintf( fp, "# FINAL MEAN DEPTH = %lf\n", (double)final_total_cov / final_nsites );
        fclose(fp);
    }
    pthread_exit(NULL);
}


/* Prints help message */
void print_help(int code) {
    fprintf( stderr, "IBDGem: Compares low-coverage sequencing data from an unknown sample to known genotypes\n" );
    fprintf( stderr, "        from a reference individual/panel and calculates the likelihood that the samples\n" );
    fprintf( stderr, "        share 0, 1, or 2 IBD chromosomes.\n\n" );
    fprintf( stderr, "Usage: ./ibdgem -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]\n" );
    fprintf( stderr, "-H, --hap  FILE                 HAP file (required)\n" );
    fprintf( stderr, "-L, --legend  FILE              LEGEND file (required)\n" );
    fprintf( stderr, "-I, --indv  FILE                INDV file (required)\n" );
    fprintf( stderr, "-P, --pileup  FILE              PILEUP file (required)\n" );
    fprintf( stderr, "-N, --pileup-name  STR          Name of individual in Pileup (default: UNKWN)\n" );
    fprintf( stderr, "-A, --allele-freqs  FILE        File containing allele frequencies from an external panel;\n" );
    fprintf( stderr, "                                   must be sorted & whitespace-delimited with columns CHROM, POS, AF;\n" );
    fprintf( stderr, "                                   use in conjunction with -c if includes multiple chromosomes\n" );
    fprintf( stderr, "                                   (default: calculate AF from input genotypes)\n" );
    fprintf( stderr, "-S, --sample-list  FILE         File containing subset of samples to compare the\n" );
    fprintf( stderr, "                                   sequencing data against; one line per sample\n" );
    fprintf( stderr, "                                   (default: compare against all samples in genotype panel)\n" );
    fprintf( stderr, "-s, --sample  STR               Sample(s) to compare the sequencing data against; comma-separated\n" );
    fprintf( stderr, "                                   without spaces if more than one (e.g. sample1,sample2,etc.)\n" );
    fprintf( stderr, "-p, --positions  FILE           List of sites to compare; can be in position list format with 2 columns\n" );
    fprintf( stderr, "                                   CHROM, POS (1-based coordinates) or BED format (0-based coordinates);\n" );
    fprintf( stderr, "                                   use in conjunction with -c if includes multiple chromosomes\n" );
    fprintf( stderr, "                                   (default: perform comparison at all sites)\n" );
    fprintf( stderr, "-M, --max-cov  INT              Maximum estimated coverage of Pileup data (default: 20)\n" );
    fprintf( stderr, "-F, --max-af  FLOAT             Maximum alternate allele frequency (default: 1)\n" );
    fprintf( stderr, "-f, --min-af  FLOAT             Minimum alternate allele frequency (default: 0)\n" );
    fprintf( stderr, "-D, --downsample-cov  FLOAT     Down-sample to this fold-coverage depth\n" );
    fprintf( stderr, "-O, --out-dir  STR              Path to output directory (default: output to current directory)\n" );
    fprintf( stderr, "-c, --chromosome  STR           Chromosome on which the comparison is done; if not specified,\n" );
    fprintf( stderr, "                                   will assume that all inputs are on one single chromosome\n" );
    fprintf( stderr, "-t, --threads  INT              Number of threads; only recommended when number of samples to compare\n" );
    fprintf( stderr, "                                   exceeds 10,000 (default: 1)\n" );
    fprintf( stderr, "-e, --error-rate  FLOAT         Error rate of sequencing platform (default: 0.02)\n" );
    fprintf( stderr, "-v, --variable-sites-only       If set, make output only for sites that have >=1\n" );
    fprintf( stderr, "                                   alternate alleles in genotype file for this sample\n" );
    fprintf( stderr, "-h, --help                      Show this help message and exit\n\n" );
    fprintf( stderr, "Format of output table is tab-delimited with columns:\n" );
    fprintf( stderr, "POS, REF, ALT, rsID, AF, DP, GT_A0, GT_A1, SQ_NREF, SQ_NALT, LIBD0, LIBD1, LIBD2\n" );
    exit(code);
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {

    clock_t start, end;
    double time_elapsed;
    start = clock();

    int option;
    char* opts = ":H:L:I:P:N:A:S:s:p:M:F:f:D:O:c:t:e:vh";

    char hap_fn[MAX_FN_LEN];
    char legend_fn[MAX_FN_LEN];
    char indv_fn[MAX_FN_LEN];
    char pu_fn[MAX_FN_LEN];
    char sample_fn[MAX_FN_LEN];
    char af_fn[MAX_FN_LEN];
    char pos_fn[MAX_FN_LEN];
    char str_cmd[MAX_FIELD_WIDTH];
    char* sample_str;
    char* out_dir;
    char* uchr = NULL;
    double target_dp = 0;
    
    if (argc == 1) {
        print_help(0);
    }
    
    char cwd[PATH_MAX];
    if ( getcwd(cwd, sizeof(cwd)) ) {
        out_dir = cwd;
    }
    else {
        // manually set output dir to current working dir
        out_dir = "./";
    }
    while ( (option = getopt_long(argc, argv, opts, longopts, NULL)) != -1 ) {
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
                SQ_ID = optarg;
                break;
            case 'S':
                OPT_S1 = 1;
                strcpy(sample_fn, optarg);
                break;
            case 's':
                OPT_S2 = 1;
                sample_str = optarg;
                break;
            case 'A':
                OPT_A = 1;
                strcpy(af_fn, optarg);
                break;
            case 'M':
                USER_MAX_COV = atoi(optarg);
                break;
            case 'F':
                USER_MAX_AF = atof(optarg);
                break;
            case 'f':
                USER_MIN_AF = atof(optarg);
                break;
            case 'p':
                OPT_P = 1;
                strcpy(pos_fn, optarg);
                break;
            case 'c':
                uchr = optarg;
                break;
            case 't':
                NUM_THREADS = atoi(optarg);
                break;
            case 'e':
                EPSILON = atof(optarg);
                break;
            case 'O':
                out_dir = optarg;
                break;
            case 'D':
                target_dp = atof(optarg);
                OPT_D = 1;
                break;
            case 'v':
                OPT_V = 1;
                break;
            case 'h':
                print_help(0);
            case ':':
                fprintf( stderr, "Option -%c missing required argument.\n", optopt );
                exit(0);
            case '?':
                if (isprint(optopt)) {
                    fprintf( stderr, "Invalid option -%c.\n", optopt );
                }
                else {
                    fprintf( stderr, "Invalid option character.\n" );
                }
                break;
            default:
                fprintf( stderr, "[::] ERROR parsing command-line options.\n" );
                exit(0);
        }
    }
    for (int i = optind; i < argc; i++) {
        fprintf( stderr, "Given extra argument %s.\n", argv[i] );
    }

    // first copy the first argument into str_cmd to initialize it
    strcpy(str_cmd, argv[0]);
    strcat(str_cmd, " ");
    for (int i = 1; i < argc; i++) {
        strcat(str_cmd, argv[i]);
        strcat(str_cmd, " ");
    }

    Impute2* i2 = init_I2(hap_fn, legend_fn, indv_fn);
    Pu_chr* puc = init_Pu_chr(pu_fn, uchr);
    if (!i2 || !puc) {
        fprintf( stderr, "[::] ERROR parsing genotype and/or sequence data; make sure input is valid.\n" );
        exit(1);
    }
    
    int n_samples = i2->n_haps/2;
    // if user provides file with names of samples to compare,
    //    overwrite i2->samples with that subset
    // if not, will do comparison for every sample in panel
    if (OPT_S1) { 
        n_samples = read_sf(i2, sample_fn);
        if (n_samples == -1) {
            destroy_I2(i2);
            destroy_Pu_chr(puc);
            exit(1);
        }
    }

    if (OPT_S2) {
        n_samples = read_scmd(i2, sample_str);
        if (n_samples == -1) {
            destroy_I2(i2);
            destroy_Pu_chr(puc);
            exit(1);
        }
    }

    if (OPT_A) {
        int af_status = read_af(i2, af_fn, uchr);
        if (af_status == 1) {
            destroy_I2(i2);
            destroy_Pu_chr(puc);
            exit(1);
        }
    }

    if (OPT_P) {
        int pos_status = read_pos(i2, pos_fn, uchr);
        if (pos_status == 1) {
            destroy_I2(i2);
            destroy_Pu_chr(puc);
            exit(1);
        }
    }
    
    if ( OPT_D && target_dp <= 0 ) {
        fprintf( stderr, "Invalid target depth of coverage (must be > 0).\n" );
        destroy_I2(i2);
        destroy_Pu_chr(puc);
        exit(0);
    }

    Prob_sum** ptab = malloc(i2->n_sites * sizeof(Prob_sum*));
    for (int i = 0; i < i2->n_sites; i++) {
        ptab[i] = malloc(sizeof(Prob_sum));
        ptab[i]->failed_filters = 0;
    }
    unsigned long** nCk = init_nCk(USER_MAX_COV);
    Cov_dist* idist = malloc(sizeof(Cov_dist));
    memset(idist->dist, 0, sizeof(unsigned long) * (USER_MAX_COV+1));
    
    double cull_p = do_pD_calc(nCk, ptab, i2, puc, target_dp, idist);

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
    Sampl* batches[NUM_THREADS];

    int cd_status = chdir(out_dir);
    if (cd_status) {
        fprintf( stderr, "[::] ERROR in do_compare(): Cannot access output directory %s.\n", out_dir );
        exit(1);
    }
    for (int ti = 0; ti < NUM_THREADS; ti++) {
        // last thread? if so will also process the remaining samples
        if (ti == NUM_THREADS-1) {
            thr_nsamples = batch_size+r;
            include_r = 1;
        }
        batches[ti] = malloc(thr_nsamples * sizeof(Sampl));
        make_batch(i2->samples, batches[ti], ti, batch_size, r, include_r);

        cargs[ti].str_cmd = str_cmd;
        cargs[ti].i2 = i2;
        cargs[ti].ptab = ptab;
        cargs[ti].samples = batches[ti];
        cargs[ti].n_samples = thr_nsamples;
        cargs[ti].cull_p = cull_p;
        cargs[ti].dist = idist;

        if ((rc = pthread_create(&thr[ti], NULL, do_compare, &cargs[ti]))) {
            fprintf( stderr, "[::] ERROR in pthread_create(), return code = %d.\n", rc);
            for (int i = 0; i < i2->n_sites; i++) {
                free(ptab[i]);
            }
            free(ptab);
            free(batches[ti]);
            free(idist);
            destroy_I2(i2);
            destroy_nCk(nCk, USER_MAX_COV);
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
    free(ptab);
    free(idist);
    destroy_I2(i2); 
    destroy_nCk(nCk, USER_MAX_COV);
    destroy_Pu_chr(puc);

    end = clock();
    time_elapsed = ( (double)(end-start)  / CLOCKS_PER_SEC ) / 60;
    fprintf( stderr, "Run time: %f minutes.\n", time_elapsed );
    return EXIT_SUCCESS;
}



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