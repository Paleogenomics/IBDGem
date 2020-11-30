#include <math.h>
#include <unistd.h> 
#include "file-io.h"
#include "load_i2.h"
#include "pileup.h"
#include "nchoosek.h"

#define EPSILON 0.02
#define BASES "ACGT"

unsigned long** N_CHOOSE_K;
double MEAN_DEPTH;

typedef struct compare {
    unsigned int* pos;
    char** ref;
    char** alt;
    char** rsID;
    double* af;
    unsigned int* dp;
    unsigned short* hap1_allele;
    unsigned short* hap2_allele;
    unsigned int* ref_counts;
    unsigned int* alt_counts;
    double* IBD0;
    double* IBD1;
    double* IBD2;
    size_t nsites_included;
} Compare;

int destroy_Compare(Compare* comp) {
    size_t n = comp->nsites_included;
    if (!comp) {
        return 0;
    }
    free_cptrarr(comp->ref, n);
    free_cptrarr(comp->alt, n);
    free_cptrarr(comp->rsID, n);
    free(comp->af);
    free(comp->dp);
    free(comp->hap1_allele);
    free(comp->hap2_allele);
    free(comp->ref_counts);
    free(comp->alt_counts);
    free(comp->IBD0);
    free(comp->IBD1);
    free(comp->IBD2);
    return 0;
}

void populate_info(Compare* comp, Impute2* i2, Pu_chr* puc, size_t sample_idx, bool opt_v) {
    size_t n = i2->num_sites;
    comp->pos = malloc(n * sizeof(unsigned int));
    comp->ref = malloc(n * sizeof(char*));
    comp->alt = malloc(n * sizeof(char*));
    comp->rsID = malloc(n * sizeof(char*));
    comp->af = malloc(n * sizeof(double));
    comp->dp = malloc(n * sizeof(unsigned int));
    comp->hap1_allele = malloc(n * sizeof(unsigned short));
    comp->hap2_allele = malloc(n * sizeof(unsigned short));
    comp->ref_counts = malloc(n * sizeof(unsigned int));
    comp->alt_counts = malloc(n * sizeof(unsigned int));

    size_t count = 0;
    for (int i = 0; i < n; i++) {
        char* ref = i2->ref_alleles[i];
        char* alt = i2->alt_alleles[i];
        unsigned short A0 = *(i2->haps[i]+sample_idx);
        unsigned short A1 = *(i2->haps[i]+(sample_idx+1));
        if (opt_v && A0 == 0 && A1 == 0) {
           continue;
        }
        else if (biallelic(ref, alt)) {
            unsigned int pos = atoi(i2->pos[i]);
            comp->pos[i] = pos;
            comp->ref[i] = ref;
            comp->alt[i] = alt;
            comp->rsID[i] = i2->ids[i];
            comp->af[i] = find_f(i2, i);
            Pul* pul = fetch_Pul(puc, pos);
            comp->dp[i] = pul->cov;
            comp->hap1_allele[i] = A0; 
            comp->hap2_allele[i] = A1;
            comp->ref_counts[i] = count_allele_from_pul(pul, ref);
            comp->alt_counts[i] = count_allele_from_pul(pul, alt);
            count++;
        }
    }
    comp->nsites_included = count;
}

void populate_likelihoods(Compare* comp) {
    for (int i = 0; i < comp->nsites_included; i++) {
        unsigned int nREF = comp->ref_counts[i];
        unsigned int nALT = comp->alt_counts[i];
        unsigned short A0 = comp->hap1_allele[i];
        unsigned short A1 = comp->hap2_allele[i];
        double f = comp->af[i];

        comp->IBD0[i] = find_LD_given_f(nREF, nALT, f);
        comp->IBD1[i] = find_LD_IBD1(A0, A1, nREF, nALT, f);
        comp->IBD2[i] = find_LD_given_G(A0, A1, nREF, nALT);
    }
}

void output_depth_dist(Compare* comp, unsigned int max_depth) {
    unsigned long total_cov = 0;
    unsigned int cov_dist[max_depth+1] = {0};
    for (int i = 0; i < comp->nsites_included; i++) {
        unsigned int cov = comp->ref_counts[i] + comp->alt_counts[i];
        if (cov <= max_depth) {
            cov_dist[cov]++;
        }
    }
    //printf( "# COVERAGE NUM_SITES\n" );
    for (int j = 0; j <= max_depth; j++) {
        //printf( "# %d %d\n", $i, $cov_dist[$i] );
        total_cov += (cov_dist[j] * j);
    }
    MEAN_DEPTH = total_cov / (double)comp->nsites_included;
}


bool biallelic(const char* ref, const char* alt) {
    if ((strstr(BASES, ref) && strlen(ref) == 1) && 
    (strstr(BASES, alt) && strlen(alt) == 1)) {
        return true;
    }
    return false;
}

int cull_depth(Compare* comp, double target_depth) {
    if (target_depth > MEAN_DEPTH) {
        fprintf(stderr, "Current mean depth of coverage is lower than -D. No culling will be done.\n");
        return 1;
    }
    unsigned int new_ref_count = 0;
    unsigned int new_alt_count = 0;
    double cull_p = target_depth / MEAN_DEPTH;
    for (int i = 0; i < comp->nsites_included; i++) {
        for (int j = 0; j < comp->ref_counts[i]; j++) {
            if(rand() / (double)RAND_MAX < cull_p) {
                new_ref_count++;
            }
        }
        comp->ref_counts[i] = new_ref_count;
        for (int j = 0; j < comp->alt_counts[i]; j++) {
            if(rand() / (double)RAND_MAX < cull_p) {
                new_alt_count++;
            }
        }
        comp->alt_counts[i] = new_alt_count;
    }
    return 0;
} 

unsigned int find_sample_idx(Impute2* i2, const char* identifier) {
    for (unsigned int i = 0; i < (i2->num_haps)/2; i++) {
        if (strcmp(i2->samples[i], identifier) == 0) {
            return i*2;
        }
    }
    fprintf(stderr, "Sample %s not found in .indv file\n", identifier);
    exit(1);
}

double find_f(Impute2* i2, size_t pos_index) {
    unsigned short* allele_arr = i2->haps[pos_index];
    unsigned int nALT = 0; 
    for (int i = 0; i < i2->num_haps; i++) {
        if (*(allele_arr+i) == 1) {
            nALT++;
        }
    }
    return nALT / i2->num_haps;
}

size_t count_allele_from_pul(Pul* pul, const char* allele) {
    size_t count = 0;
    for (int i; i < pul->cov; i++) {
        if (pul->bases[i] == allele) {
            count++;
        }
    }
    return count;
}

double find_LD_given_G(unsigned short A0, unsigned short A1, size_t nREF, size_t nALT) {
    double L = 1;
    if (nREF == 0 && nALT == 0) {
        return 1;
    }
    unsigned long N_choose_k = retrieve(N_CHOOSE_K, nREF+nALT, nREF);
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
    return L;
}

double find_LD_given_f(size_t nREF, size_t nALT, double f) {
    double L = 1;
    if (nREF == 0 && nALT == 0) {
        return 1;
    }
    L = pow(1-f, 2.0) * find_LD_given_G(0, 0, nREF, nALT) +
        2 * (1-f) * f * find_LD_given_G(0, 1, nREF, nALT) +
        pow(f, 2.0) * find_LD_given_G(1, 1, nREF, nALT);
    return L;
}

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

int main(int argc, char* argv[]) {
    int option;
    char* opts = ":H:L:I:P:S:M:D:v";
    bool opt_v = false;
    bool opt_D = false;
    char* hap_fn, legend_fn, indv_fn, pu_fn, sample_id;
    unsigned int max_depth = 0;
    unsigned int target_cull_depth = 0;
    while ((option = getopt(argc, argv, opts)) != -1) {
        switch (option) {
            case 'H':
                hap_fn = optarg;
                break;
            case 'L':
                legend_fn = optarg;
                break;
            case 'I':
                indv_fn = optarg;
                break;
            case 'P':
                pu_fn = optarg;
                break;
            case 'S':
                sample_id = optarg;
                break;
            case 'M':
                max_depth = atoi(optarg);
                break;
            case 'D':
                target_cull_depth = atoi(optarg);
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
                exit(0);
            default:
                fprintf(stderr, "Error parsing command-line options.\n");
                exit(0);
        }
        for (int i = optind; i < argc; i++) {
            printf("Non-option argument %s\n", argv[i]);
        }
    }
    if (strlen(sample_id) == 0) {
        fprintf(stderr, "Please specify sample ID via the -S option.\n");
        exit(0);
    }
    if (max_depth == 0) {
        fprintf(stderr, "Please specify maximum depth of coverage via the -M option.\n");
        exit(0);
    }
    if (opt_D && target_cull_depth == 0) {
        fprintf(stderr, "Please specify culling depth via the -D option.\n");
        exit(0);
    }
    Impute2* i2 = init_I2(hap_fn, legend_fn, indv_fn);
    Pu_chr* puc = init_Pu_chr(pu_fn);
    Compare* comp = malloc(sizeof(Compare));
    N_CHOOSE_K = init_nCk(max_depth);
    unsigned int sample_idx = find_sample_idx(i2, sample_id);
    populate_info(comp, i2, puc, sample_id, opt_v);
    output_depth_dist(comp, max_depth);
    if (opt_D && target_cull_depth > 0) {
        cull_depth(comp, target_cull_depth);
        output_depth_dist(comp, max_depth);
    }
    populate_likelihoods(comp);
    //call function for writing results to file
    destroy_Compare(comp);
    destroy_I2(i2);
    destroy_nCk(N_CHOOSE_K, max_depth);
    //destroy Pu_chr
    return 0;
}

