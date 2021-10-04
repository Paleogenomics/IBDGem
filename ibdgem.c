#include <math.h>
#include <unistd.h>
#include <time.h>
#include <float.h>

#include "file-io.h"
#include "load_i2.h"
#include "pileup.h"
#include "nchoosek.h"

#define EPSILON 0.02
#define BASES "ACGT"
unsigned int USER_MAX_DEPTH = 20;
unsigned long** N_CHOOSE_K;


size_t count_allele_from_pul(Pul* pul, const char allele) {
    size_t count = 0;
    for (int i = 0; i < pul->cov; i++) {
        if (pul->bases[i] == allele) {
            count++;
        }
    }
    return count;
}

double find_cull_p(bool opt_D, Pu_chr* puc, double target_depth) {
    double cull_p;
    size_t cov_dist[USER_MAX_DEPTH+1];
    for (int i = 0; i <= USER_MAX_DEPTH; i++) {
        cov_dist[i] = 0;
    }
    size_t total_cov = 0;
    for (int j = 0; j < puc->n_puls; j++) {
        unsigned int cov = puc->puls[j]->cov;
        if (cov <= USER_MAX_DEPTH) {
            total_cov += cov;
            cov_dist[cov]++;
        }
    }
    printf("# INITAL COVERAGE DISTRIBUTION:\n");
    printf("# COVERAGE NUM_SITES\n");
    for (int n = 0; n <= USER_MAX_DEPTH; n++) {
        printf("# %d %zu\n", n, cov_dist[n]);
    }
    double mean_depth = (double)total_cov / puc->n_puls;
    if (target_depth > mean_depth) {
        fprintf(stderr, "Observed depth is lower than target depth -D. No culling will be done.\n");
        return 1;
    }
    printf("# MEAN DEPTH = %lf\n", mean_depth);
    if  (!opt_D) {
        cull_p = 1.00;
    }
    else {
        cull_p = target_depth / mean_depth;
    }
    printf("# CULL DEPTH RATIO = %lf\n", cull_p);
    return cull_p;
}

bool biallelic(const char* ref, const char* alt) {
    if ((strstr(BASES, ref) && strlen(ref) == 1) && 
    (strstr(BASES, alt) && strlen(alt) == 1)) {
        return true;
    }
    return false;
}

size_t find_sample_idx(Impute2* i2, const char* identifier) {
    for (int i = 0; i < (i2->num_haps)/2; i++) {
        if (strcmp(i2->samples[i], identifier) == 0) {
            return i*2;
        }
    }
    fprintf(stderr, "Sample %s not found in .indv file\n", identifier);
    exit(1);
}

double find_f(Impute2* i2, size_t pos_index) {
    unsigned short* allele_arr = i2->haps[pos_index];
    size_t nALT = 0; 
    for (int i = 0; i < i2->num_haps; i++) {
        if (*(allele_arr+i) == 1) {
            nALT++;
        }
    } 
    return (double)nALT / i2->num_haps;
}

size_t cull_depth(unsigned int original_count, double cull_p) {
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

double find_LD_given_G(unsigned short A0, unsigned short A1, size_t nREF, size_t nALT) {
    double L = 1.00; 
    
    if (nREF == 0 && nALT == 0) {
        return L;
    }
    unsigned long N_choose_k = retrieve_nCk(N_CHOOSE_K, nREF+nALT, nREF);
    if (A0 == 0 && A1 == 0) {
        L = N_choose_k * pow(1-EPSILON, nREF) * pow(EPSILON, nALT);
    }
    else if (A0 == 1 && A1 == 1) {
        L = N_choose_k * pow(1-EPSILON, nALT) * pow(EPSILON, nREF);
    }
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        L = N_choose_k * pow(0.5, nREF) * pow(0.5, nALT);
    }
    else {
        fprintf(stderr, "Invalid genotype: A0 = %d, A1 = %d\n", A0, A1);
        exit(1);
    }

    if (L == 0.00) {
        L = DBL_MIN;
    }
    return L;
}

double find_LD_given_f(size_t nREF, size_t nALT, double f) {
    double L = 1.00;
    
    if (nREF == 0 && nALT == 0) {
        return L;
    }
    L = pow(1-f, 2.0) * find_LD_given_G(0, 0, nREF, nALT) +
        2 * (1-f) * f * find_LD_given_G(0, 1, nREF, nALT) +
        pow(f, 2.0) * find_LD_given_G(1, 1, nREF, nALT);

    if (L == 0.00) {
        L = DBL_MIN;
    }
    return L;
}

double find_LD_IBD1(unsigned short A0, unsigned short A1, size_t nREF, size_t nALT, double f) {
    double L = 1.00;

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

    if (L == 0.00) {
        L = DBL_MIN;
    }
    return L;
}
void output_table(Impute2* i2, Pu_chr* puc, const char* sample_id, double target_depth, bool opt_D, bool opt_v) {
    size_t sample_idx = find_sample_idx(i2, sample_id);
    double cull_p = find_cull_p(opt_D, puc, target_depth);
    
    size_t final_nsites = 0;
    size_t final_total_cov = 0;
    size_t final_cov_dist[USER_MAX_DEPTH+1];
    for (int cov = 0; cov <= USER_MAX_DEPTH; cov++) {
        final_cov_dist[cov] = 0;
    }
    printf("# Pos\tREF\tALT\trsID\tAF\tDP\tVCFA0\tVCFA1\tBAMnREF\tBAMnALT\tL(IBD0)\tL(IBD1)\tL(IBD2)\n");
    for (int i = 0; i < i2->num_sites; i++) {
        char* ref = i2->ref_alleles[i];
        char* alt = i2->alt_alleles[i];
        unsigned short A0 = *(i2->haps[i]+sample_idx);
        unsigned short A1 = *(i2->haps[i]+(sample_idx+1));
        if (opt_v && A0 == 0 && A1 == 0) {
           continue;
        }
        else if (biallelic(ref, alt)) {
            size_t pos = atoi(i2->pos[i]);
            char* rsID = i2->ids[i];
            double f = find_f(i2, i);
            Pul* pul = fetch_Pul(puc, pos);
            if (pul) {
                unsigned int dp = pul->cov;
                size_t nREF = count_allele_from_pul(pul, ref[0]);
                size_t nALT = count_allele_from_pul(pul, alt[0]);
	    	    if (nREF + nALT <= USER_MAX_DEPTH) {
                    nREF = cull_depth(nREF, cull_p);
                    nALT = cull_depth(nALT, cull_p);
                    double IBD0 = find_LD_given_f(nREF, nALT, f);
                    double IBD1 = find_LD_IBD1(A0, A1, nREF, nALT, f);
                    double IBD2 = find_LD_given_G(A0, A1, nREF, nALT);
                    printf("%zu\t%s\t%s\t%s\t%lf\t%u\t%u\t%u\t%zu\t%zu\t%e\t%e\t%e\n",
                        pos, ref, alt, rsID, f, dp, A0, A1, nREF, nALT, IBD0, IBD1, IBD2);
                    final_cov_dist[nREF+nALT]++;
                    final_total_cov += (nREF+nALT);
		            final_nsites++;
                }
            }
        }
    }
    printf("# FINAL COVERAGE DISTRIBUTION:\n");
    printf("# COVERAGE NUM_SITES\n");
    for (int n = 0; n <= USER_MAX_DEPTH; n++) {
        printf("# %d %zu\n", n, final_cov_dist[n]);
    }
    printf("# FINAL MEAN DEPTH = %lf\n", (double)final_total_cov / final_nsites);
}

int main(int argc, char* argv[]) {
    srand(time(NULL));
    int option;
    char* opts = ":H:L:I:P:S:M:D:v";
    bool opt_v = false;
    bool opt_D = false;
    char* hap_fn, *legend_fn, *indv_fn, *pu_fn, *sample_id = NULL;
    double target_cull_depth = 0;
    while ((option = getopt(argc, argv, opts)) != -1) {
        switch (option) {
            case 'H':
                hap_fn = strdup(optarg);
                break;
            case 'L':
                legend_fn = strdup(optarg);
                break;
            case 'I':
                indv_fn = strdup(optarg);
                break;
            case 'P':
                pu_fn = strdup(optarg);
                break;
            case 'S':
                sample_id = strdup(optarg);
                break;
            case 'M':
                USER_MAX_DEPTH = atoi(optarg);
                break;
            case 'D':
                target_cull_depth = atof(optarg);
                opt_D = true;
                break;
            case 'v':
                opt_v = true;
                break;
            case ':':
                fprintf(stderr, "Please enter required argument for option -%c.\n", optopt);
                exit(0);
            case '?':
                if (isprint(optopt)) {
                    fprintf(stderr, "Unknown option -%c.\n", optopt);
                }
                else {
                    fprintf (stderr, "Unknown option character \\x%x.\n", optopt);
                }
                break;
            default:
                fprintf(stderr, "Error parsing command-line options.\n");
                exit(0);
        }
    }
    for (int i = optind; i < argc; i++) {
            printf("Non-option argument %s\n", argv[i]);
    }
    if (!hap_fn || !legend_fn || !indv_fn) {
        fprintf(stderr, "Find 0, 1, or 2 IBD segments between\n" );
	    fprintf(stderr, "IMPUTE2 genotype files and BAM file info.\n" );
	    fprintf(stderr, "-H <HAP file>\n");
        fprintf(stderr, "-L <LEGEND file>\n");
        fprintf(stderr, "-I <IDNV file>\n");
	    fprintf(stderr, "-P <MPILEUP file>\n" );
	    fprintf(stderr, "-D <downsample to this fold-coverage depth>\n");
	    fprintf(stderr, "-S <identifier of sample if there are multiple samples in the INDV file>\n");
	    fprintf(stderr, "-v <if set, make output only for sites that have >=1 variant\n");
	    fprintf(stderr, "    allele in genotype file for this sample>\n");
	    fprintf(stderr, "Format of output table is tab-delimited with columns:\n");
	    fprintf(stderr, "Position, REF_ALLELE, ALT_ALLELE, rsID, AF, DP, VCFA0, VCFA1, BAMnREF, BAMnALT, L(IBD0), L(IBD1), L(IBD2)\n");
        exit(0);
    }
    if (!sample_id) {
        fprintf(stderr, "Please specify sample ID via the -S option.\n");
        exit(0);
    }
    if (USER_MAX_DEPTH == 0) {
        fprintf(stderr, "Please specify maximum depth of coverage via the -M option.\n");
        exit(0);
    }
    if (opt_D && target_cull_depth == 0) {
        fprintf(stderr, "Please specify culling depth via the -D option.\n");
        exit(0);
    }
    
    Impute2* i2 = init_I2(hap_fn, legend_fn, indv_fn);
    Pu_chr* puc = init_Pu_chr(pu_fn);
    N_CHOOSE_K = init_nCk(USER_MAX_DEPTH);
    
    output_table(i2, puc, sample_id, target_cull_depth, opt_D, opt_v);
    
    free(hap_fn);
    free(legend_fn);
    free(indv_fn);
    free(pu_fn);
    free(sample_id);
    destroy_I2(i2);
    destroy_nCk(N_CHOOSE_K, USER_MAX_DEPTH);
    destroy_Pu_chr(puc);
    return 0;
}
