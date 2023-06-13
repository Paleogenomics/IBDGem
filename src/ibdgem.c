/*** 
 * IBDGEM - Program for genetic identification from low-coverage sequencing data
 * (c) UC Regents
 * Compares sequence alignment information from an unknown sample to known genotypes
 * generated through deep-sequencing or other means and evaluates the likelihood that
 * the sample in the sequencing data shares 0, 1, or 2 IBD chromosomes with the
 * individual(s) providing the genotype data.
 *  ***/

#include <unistd.h>
#include <getopt.h>
#include <time.h>

#include "file-io.h"
#include "ibd-parse.h"
#include "pileup.h"
#include "ibd-math.h"

#define BASES "ACGT"

static double EPSILON = 0.02;
static double USER_MIN_QUAL = 0;
static unsigned int USER_MAX_COV = 20;
static double USER_MAX_AF = 1;
static double USER_MIN_AF = 0;
static int WINDOWSIZE = 100;
static char* SQ_ID = "UNKWN";

static int LD_MODE = 0;
static int IN_IMPUTE = 0; // flag for IMPUTE input
static int IN_VCF = 0; // flag for VCF input
static int OPT_S1 = 0; // flag for comparing a subset of samples from user-specified file
static int OPT_S2 = 0; // flag for comparing a subset of samples supplied through command line
static int OPT_B = 0; // flag for using user-specified samples as background panel
static int OPT_A = 0; // flag for using user-specified allele frequencies
static int OPT_P = 0; // flag for using user-specified sites
static int OPT_V = 0; // flag for skipping sites that are homozygous REF
static int OPT_D = 0; // down-sampling flag

/** Long options table **/
static struct option longopts[] = {
    { "LD",                   no_argument, &LD_MODE, 1},
    { "vcf",                  required_argument, 0, 'V' },
    { "hap",                  required_argument, 0, 'H' },
    { "legend",               required_argument, 0, 'L' },
    { "indv",                 required_argument, 0, 'I' },
    { "pileup",               required_argument, 0, 'P' },
    { "pileup-name",          required_argument, 0, 'N' },
    { "window-size",          required_argument, 0, 'w' },
    { "allele-freqs",         required_argument, 0, 'A' },
    { "sample-list",          required_argument, 0, 'S' },
    { "sample",               required_argument, 0, 's' },
    { "background-list",      required_argument, 0, 'B' },
    { "max-cov",              required_argument, 0, 'M' },
    { "downsample-cov",       required_argument, 0, 'D' },
    { "out-dir",              required_argument, 0, 'O' },
    { "max-af",               required_argument, 0, 'F' },
    { "min-af",               required_argument, 0, 'f' },
    { "positions",            required_argument, 0, 'p' },
    { "min-qual",             required_argument, 0, 'q' },
    { "chromosome",           required_argument, 0, 'c' },
    { "error-rate",           required_argument, 0, 'e' },
    { "variable-sites-only",  no_argument, 0, 'v' },
    { "help",                 no_argument, 0, 'h' },
    { 0, 0, 0, 0}
};



/* Structure for storing coverage distribution info */
typedef struct dpdist {
    unsigned long dist[MAX_COV];
    double mean_cov;
} Dpdist;


/* Calculates ratio for down-sampling sequence data to the target 
   depth of coverage
   Args: Pu_chr* puc - pointer to Pileup data
         double target_dp - coverage after down-sampling
         Dpdist* idist - input coverage distribution
   Returns: ratio for down-sampling */
double find_cull_p(Pu_chr* puc, double target_dp, Dpdist* pu_dist) {
    double cull_p = 1;
    unsigned long total_cov = 0;

    // update input distribution with coverage from Pileup line
    for (int i = 0; i < puc->n_puls; i++) {
        unsigned int cov = puc->puls[i]->cov;
        if (cov <= USER_MAX_COV) {
            total_cov += cov;
            pu_dist->dist[cov]++;
        }
    }
    double mean_cov = (double)total_cov / puc->n_puls;
    pu_dist->mean_cov = mean_cov;
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
   Args: const char* ref - reference allele at this site 
         const char* alt - alternate allele at this site
   Returns: 1 if site is single-polymorphic (not indel); 0 otherwise */
int is_snp(const char* ref, const char* alt) {
    if ((strstr(BASES, ref) && strlen(ref) == 1) && 
        (strstr(BASES, alt) && strlen(alt) == 1)) {
        return 1;
    }
    return 0;
}


/* Performs down-sampling of alignment data at a single site
   Args: unsigned int original_count - initial base coverage at site 
         double cull_p - ratio for down-sampling
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


/* Prints input coverage distribution of the sequence data
   Args: FILE* out - pointer to output file
         Dpdist* pu_dist - input coverage distribution
         double cull_p - ratio for down-sampling */
void print_pu_dist(FILE* out, Dpdist* pu_dist, double cull_p) {
    fprintf( out, "# INPUT COVERAGE DISTRIBUTION:\n" );
    fprintf( out, "# COVERAGE N_SITES\n" );
    for (int cov = 0; cov <= USER_MAX_COV; cov++) {
        fprintf( out, "# %d %lu\n", cov, pu_dist->dist[cov] );
    }
    fprintf( out, "# MEAN DEPTH = %lf\n", pu_dist->mean_cov );
    fprintf( out, "# CULL DEPTH RATIO = %lf\n", cull_p );
}


/* Checks if site is biallelic (VCF)
   Args: const char* alt - pointer to alternate allele string 
   Returns: 1 if biallelic, 0 otherwise */
int is_biallelic(const char* alt) {
    // see if alt contains a comma
    if (strchr(alt, ',') == NULL) {
        return 1;
    }
    return 0;
}


/* Performs comparison using VCF input and outputs the results to 2 files:
   (1) Table file (*.tab.txt) with information about each site and their
       associated IBD0, IBD1, and IBD2 likelihoods (these likelihoods
       are NOT calculated under LD mode, even with --LD option specified)
   (2) Summary file (*.summary.txt) with likelihoods for models IBD0, IBD1,
       and IBD2 aggregated over multiple sites, determined by --window-size
       (these aggregated likelihoods would be calculated under LD mode
       if --LD option was specified)

   Args: File_Src* vcf_fp - pointer to VCF file
         Comp_dt* data - pointer to comparison data
         Pu_chr* puc - pointer to Pileup data
         unsigned long** nCk - pointer to the nchoosek matrix
         double target_dp - coverage after down-sampling
         Dpdist* pu_dist - input coverage distribution
         const char* out_dir - path to output directory
         const char* user_cmd - entered user command
   Returns: 0 on success; 1 if error */
int compare_vcf(File_Src* vcf_fp, Comp_dt* data, Pu_chr* puc, unsigned long** nCk,
                double target_dp, Dpdist* pu_dist, const char* out_dir, const char* user_cmd) {

    int pu_idx = -1;
    Sampl* pu_sample = find_sample(data, SQ_ID);
    if (pu_sample) {
        pu_idx = pu_sample->idx;
    }

    double cull_p = find_cull_p(puc, target_dp, pu_dist);
    if (!data->uids) {
        data->uids = data->ids;
        data->n_uids = data->n_ids;
    }
    if (!data->refids) {
        data->refids = data->ids;
        data->n_refids = data->n_ids;
    }

    regex_t regex;
    int rc_flag = regcomp(&regex, "^[01][/|][01].*$", 0);
    if (rc_flag) {
        fprintf( stderr, "[::] ERROR in read_vcf(): Cannot compile regex to check genotype.\n" );
        return 1;
    }

    for (int i = 0; i < data->n_uids; i++) {
        fprintf( stderr, "Running %s-vs-%s comparison...\n", SQ_ID, data->uids[i].name );        
        char line[MAX_LINE_LEN];
        get_line_FS(vcf_fp, line);
        while (strstr(line, "##") == line) {
            get_line_FS(vcf_fp, line);
        }

        char tab_fn[MAX_FN_LEN];
        char sum_fn[MAX_FN_LEN];
        sprintf( tab_fn, "%s/%s.%s.tab.txt", out_dir, SQ_ID, data->uids[i].name );
        sprintf( sum_fn, "%s/%s.%s.summary.txt", out_dir, SQ_ID, data->uids[i].name );
        FILE* tab_fp = fopen(tab_fn, "w");
        FILE* sum_fp = fopen(sum_fn, "w");
        if (!tab_fp || !sum_fp) {
            fprintf( stderr, "[::] ERROR in compare_vcf(): Cannot open '%s' and/or '%s' for writing.\n", tab_fn, sum_fn );
            return 1;
        }
        fprintf( tab_fp, "# Entered command: %s\n\n", user_cmd );
        
        unsigned long processed = 0;
        unsigned long skipped = 0;
        print_pu_dist(tab_fp, pu_dist, cull_p);
        unsigned long final_dist[USER_MAX_COV+1];
        for (int cov = 0; cov < USER_MAX_COV+1; cov++) {
            final_dist[cov] = 0;
        }
        unsigned long final_total_cov = 0;

        fprintf( tab_fp, "# CHR\trsID\tPOS\tREF\tALT\tAF\tDP\tSQ_NREF\tSQ_NALT\tGT_A0\tGT_A1\tLIBD0\tLIBD1\tLIBD2\n" );
        fprintf( sum_fp, "# SEGMENT\tSTART\tEND\tLIBD0\tLIBD1\tLIBD2\tNUM_SITES\n" );

        unsigned int cmp_idx = data->uids[i].idx;

        int sgmt_count = 1;
        get_line_FS(vcf_fp, line);
        char* read_vcf_res = &line[0];
        while (read_vcf_res) {
            int snp_count = 0;
            unsigned long sgmt_start, sgmt_end, previous_pos;

            double sum_ibd0 = 1, sum_ibd1 = 1, sum_ibd2 = 1;
            double sum_ibd2_ref[data->n_refids];
            for (int j = 0; j < data->n_refids; j++) {
                sum_ibd2_ref[j] = 1;
            }

            while (snp_count < WINDOWSIZE) {
                read_vcf_res = get_line_FS(vcf_fp, line);
                if (!read_vcf_res) {
                    sgmt_end = previous_pos;
                    break;
                }
                unsigned long pos;
                double f;
                char id[512], ref[512], alt[512], qual[512], fltr[512], info[7680];
                char gt[MAX_LINE_LEN];
                
                if (sscanf(line, "%*s %lu %512[^\t] %512[^\t] %512[^\t] %512[^\t] %512[^\t] %7680[^\t] %30720[^\n]",
                    &pos, id, ref, alt, qual, fltr, info, gt) == 8) {
                    // skip site if not biallelic        
                    if (!is_biallelic(alt)) {
                        skipped++;
                        continue;
                    }
                    unsigned short* alleles = malloc((data->n_ids*2) * sizeof(unsigned short));
                    int parse_gt_res = vcf_parse_gt(gt, regex, data->n_ids, alleles);
                    if (parse_gt_res) {
                        fprintf( stderr, "Failed to parse genotype fields at %lu. Skipping to next site.\n", pos );
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    if (OPT_V && alleles[cmp_idx] == 0 && alleles[cmp_idx+1] == 0) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    if ( !is_snp(ref, alt) ) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    if ( atof(qual) < USER_MIN_QUAL ) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    Pul* pul = fetch_Pul(puc, pos);
                    if (!pul) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    if (OPT_P && !site_in_upos(data, pos)) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    f = find_f_vcf(alleles, data->n_ids);
                    if (OPT_A) {
                        Freq* fp = fetch_freq(data, pos);
                        if (fp) {
                            f = fp->f;
                        }
                    }
                    if (f > USER_MAX_AF || f < USER_MIN_AF) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    unsigned int n_ref = count_base_from_pul(pul, ref[0]);
                    unsigned int n_alt = count_base_from_pul(pul, alt[0]);
                    if (n_ref+n_alt > USER_MAX_COV) {
                        free(alleles);
                        skipped++;
                        continue;
                    }
                    n_ref = cull_dp(n_ref, cull_p);
                    n_alt = cull_dp(n_alt, cull_p);
                    final_total_cov += (n_ref + n_alt);
                    final_dist[n_ref+n_alt]++;

                    double pDg00 = find_pDgG(nCk, EPSILON, 0, 0, n_ref, n_alt);
                    double pDg01 = find_pDgG(nCk, EPSILON, 0, 1, n_ref, n_alt);
                    double pDg11 = find_pDgG(nCk, EPSILON, 1, 1, n_ref, n_alt);
                
                    unsigned short A0, A1;
                    double ibd0 = 1, ibd1 = 1, ibd2 = 1;
                    A0 = alleles[cmp_idx];
                    A1 = alleles[cmp_idx+1];
                
                    ibd0 = find_pDgf(f, pDg00, pDg01, pDg11);
                    ibd1 = find_pDgIBD1(A0, A1, f, pDg00, pDg01, pDg11);
                    if (A0 == 0 && A1 == 0) {
                        ibd2 = pDg00;
                    }
                    else if ( (A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0) ) {
                        ibd2 = pDg01;
                    }
                    else if (A0 == 1 && A1 == 1) {
                        ibd2 = pDg11;
                    }
                    if (n_ref+n_alt < 1) {
                        fprintf( tab_fp, "%s\t%s\t%lu\t%s\t%s\t%lf\t%u\t%u\t%u\t%hu\t%hu\t%e\t%e\t%e\n", 
                                 pul->chr, id, pos, ref, alt, f, pul->cov, n_ref, n_alt, 
                                 A0, A1, ibd0, ibd1, ibd2 ); 
                        processed++;
                        free(alleles);
                        continue;
                    }
                
                    sum_ibd0 *= ibd0;
                    sum_ibd1 *= ibd1;
                    sum_ibd2 *= ibd2;

                    if (LD_MODE) {
                        for (int n = 0; n < data->n_refids; n++) {
                            double ref_ibd2;
                            int ref_idx = data->refids[n].idx;
                            A0 = alleles[n*2];
                            A1 = alleles[(n*2)+1];
                            if (A0 == 0 && A1 == 0) {
                                ref_ibd2 = pDg00;
                            }
                            else if ( (A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0) ) {
                                ref_ibd2 = pDg01;
                            }
                            else if (A0 == 1 && A1 == 1) {
                                ref_ibd2 = pDg11;
                            }
                            if (ref_idx != pu_idx) {
                                sum_ibd2_ref[n] *= ref_ibd2;
                            }
                        }
                    }
                    snp_count++;
                    previous_pos = pos;
                    if (snp_count == 1) {
                        sgmt_start = pos;
                    }
                    else if (snp_count == WINDOWSIZE) {
                        sgmt_end = previous_pos;
                    }
                    fprintf( tab_fp, "%s\t%s\t%lu\t%s\t%s\t%lf\t%u\t%u\t%u\t%hu\t%hu\t%e\t%e\t%e\n", 
                             pul->chr, id, pos, ref, alt, f, pul->cov, n_ref, n_alt, 
                             alleles[cmp_idx], alleles[cmp_idx+1], ibd0, ibd1, ibd2 ); 
                    processed++;
                    free(alleles);
                }   
                else {
                    skipped++;
                }
            }
            int n_refpanel = data->n_refids;
            if (LD_MODE) {
                double sum_ibd0_LD = 0;
                for (int k = 0; k < data->n_refids; k++) {
                    if (data->refids[k].idx != pu_idx) {
                        sum_ibd0_LD += sum_ibd2_ref[k];
                    }
                    else {
                        n_refpanel--;
                    }
                }
                fprintf( sum_fp, "%d\t%lu\t%lu\t%e\t%e\t%e\t%d\n", sgmt_count, sgmt_start,
                         sgmt_end, sum_ibd0_LD/n_refpanel, sum_ibd1, sum_ibd2, snp_count );
            }
            else {
                fprintf( sum_fp, "%d\t%lu\t%lu\t%e\t%e\t%e\t%d\n", sgmt_count, sgmt_start,
                         sgmt_end, sum_ibd0, sum_ibd1, sum_ibd2, snp_count );
            }
            sgmt_count++;
        }    
        fprintf( tab_fp, "# FINAL COVERAGE DISTRIBUTION:\n" );
        fprintf( tab_fp, "# COVERAGE N_SITES\n" );
        for (int cov = 0; cov < USER_MAX_COV+1; cov++) {
            fprintf( tab_fp, "# %d %lu\n", cov, final_dist[cov] );
        }
        fprintf( tab_fp, "# FINAL MEAN DEPTH = %lf\n", (double)final_total_cov/processed );
        fprintf( tab_fp , "## Number of sites processed: %lu\n", processed);
        fprintf( tab_fp, "## Number of sites skipped: %lu\n", skipped);
        fclose(tab_fp);
        fclose(sum_fp);
        regfree(&regex);
        rewind_FS(vcf_fp);
    }
    return 0;
}
            

/* Performs comparison using IMPUTE input and outputs the results to 2 files:
   (1) Table file (*.tab.txt) with information about each site and their
       associated IBD0, IBD1, and IBD2 likelihoods (these likelihoods
       are NOT calculated under LD mode, even with --LD option specified)
   (2) Summary file (*.summary.txt) with likelihoods for models IBD0, IBD1,
       and IBD2 aggregated over multiple sites, determined by --window-size
       (these aggregated likelihoods would be calculated under LD mode
       if --LD option was specified)

   Args: File_Src* hap_fp - pointer to .hap file
         File_Src* legend_fp - pointer to .legend file
         Comp_dt* data - pointer to comparison data
         Pu_chr* puc - pointer to Pileup data
         unsigned long** nCk - pointer to the nchoosek matrix
         double target_dp - coverage after down-sampling
         Dpdist* pu_dist - input coverage distribution
         const char* out_dir - path to output directory
         const char* user_cmd - entered user command
   Returns: 0 on success; 1 if error */
int compare_impute(File_Src* hap_fp, File_Src* legend_fp, Comp_dt* data, Pu_chr* puc, unsigned long** nCk,
                   double target_dp, Dpdist* pu_dist, const char* out_dir, const char* user_cmd) {

    int pu_idx = -1;
    // check to see if the subject individual is also in reference panel
    Sampl* pu_sample = find_sample(data, SQ_ID);
    if (pu_sample) {
        pu_idx = pu_sample->idx;
    }

    double cull_p = find_cull_p(puc, target_dp, pu_dist);
    // if subset of samples to compare not specified,
    // compare against all individuals in VCF
    if (!data->uids) {
        data->uids = data->ids;
        data->n_uids = data->n_ids;
    }
    // if reference panel not specified, use
    // all individuals in VCF as reference
    if (!data->refids) {
        data->refids = data->ids;
        data->n_refids = data->n_ids;
    }

    for (int i = 0; i < data->n_uids; i++) {
        fprintf( stderr, "Running %s-vs-%s comparison...\n", SQ_ID, data->uids[i].name );
        
        char tab_fn[MAX_FN_LEN]; // output file with IBD likelihoods per site
        char sum_fn[MAX_FN_LEN];  // output file with IBD likelihoods per genomic segment
        sprintf( tab_fn, "%s/%s.%s.tab.txt", out_dir, SQ_ID, data->uids[i].name );
        sprintf( sum_fn, "%s/%s.%s.summary.txt", out_dir, SQ_ID, data->uids[i].name );
        FILE* tab_fp = fopen(tab_fn, "w");
        FILE* sum_fp = fopen(sum_fn, "w");
        if (!tab_fp || !sum_fp) {
            fprintf( stderr, "[::] ERROR in compare_impute(): Cannot open '%s' and/or '%s' for writing.\n", tab_fn, sum_fn );
            return 1;
        }
        fprintf( tab_fp, "# Entered command: %s\n\n", user_cmd );
        
        unsigned long processed = 0;
        unsigned long skipped = 0;
        // print original coverage distribution of sequence data
        print_pu_dist(tab_fp, pu_dist, cull_p);
        unsigned long final_dist[USER_MAX_COV+1];
        for (int cov = 0; cov < USER_MAX_COV+1; cov++) {
            final_dist[cov] = 0;
        }
        unsigned long final_total_cov = 0;
        
        fprintf( tab_fp, "# CHR\trsID\tPOS\tREF\tALT\tAF\tDP\tSQ_NREF\tSQ_NALT\tGT_A0\tGT_A1\tLIBD0\tLIBD1\tLIBD2\n" );
        fprintf( sum_fp, "# SEGMENT\tSTART\tEND\tLIBD0\tLIBD1\tLIBD2\tNUM_SITES\n" );

        unsigned int cmp_idx = data->uids[i].idx; // index position of compared sample
        char hap_buf[MAX_LINE_LEN];
        char legend_buf[MAX_LINE_LEN];
        
        int sgmt_count = 1;
        get_line_FS(legend_fp, legend_buf); // skip header in legend file
        char* read_hap_res = &hap_buf[0];
        char* read_legend_res = &legend_buf[0];
        while (read_hap_res && read_legend_res) {
            int snp_count = 0;
            unsigned long sgmt_start, sgmt_end, previous_pos;
            
            double sum_ibd0 = 1, sum_ibd1 = 1, sum_ibd2 = 1;
            double sum_ibd2_ref[data->n_refids];
            for (int j = 0; j < data->n_refids; j++) {
                sum_ibd2_ref[j] = 1;
            }

            while (snp_count < WINDOWSIZE) {
                read_hap_res = get_line_FS(hap_fp, hap_buf);
                read_legend_res = get_line_FS(legend_fp, legend_buf);
                if (!read_hap_res || !read_legend_res) {
                    sgmt_end = previous_pos;
                    break;
                }
                unsigned long pos;
                double f;
                char id[128], ref[128], alt[128];
                
                // check if homozygous ref in compared sample
                if (OPT_V && hap_buf[cmp_idx*2] == '0' && hap_buf[(cmp_idx*2)+2] == '0') {
                    skipped++;
                    continue;
                }
                // check if legend line looks good
                if ( sscanf(legend_buf, "%128s %lu %128s %128s", id, &pos, ref, alt) != 4 ) {
                    skipped++;
                    continue;
                }
                if ( !is_snp(ref, alt) ) { // check if site is SNP
                    skipped++;
                    continue;
                }
                // check if site has data in Pileup
                Pul* pul = fetch_Pul(puc, pos);
                if (!pul) {
                    skipped++;
                    continue;
                }
                // check if site is in user-specified list
                if (OPT_P && !site_in_upos(data, pos)) {
                    skipped++;
                    continue;
                }
                f = find_f_impute(hap_buf, data->n_ids);
                if (OPT_A) {
                    Freq* fp = fetch_freq(data, pos);
                    if (fp) {
                        f = fp->f;
                    }
                }
                // check if AF is within specified range
                if (f > USER_MAX_AF || f < USER_MIN_AF) {
                    skipped++;
                    continue;
                }
                unsigned int n_ref = count_base_from_pul(pul, ref[0]);
                unsigned int n_alt = count_base_from_pul(pul, alt[0]);
                // check if site coverage is under specified maximum
                if (n_ref+n_alt > USER_MAX_COV) {
                    skipped++;
                    continue;
                }
                n_ref = cull_dp(n_ref, cull_p);
                n_alt = cull_dp(n_alt, cull_p);
                final_total_cov += (n_ref + n_alt);
                final_dist[n_ref+n_alt]++;

                double pDg00 = find_pDgG(nCk, EPSILON, 0, 0, n_ref, n_alt);
                double pDg01 = find_pDgG(nCk, EPSILON, 0, 1, n_ref, n_alt);
                double pDg11 = find_pDgG(nCk, EPSILON, 1, 1, n_ref, n_alt);
                
                unsigned short A0, A1;
                double ibd0 = 1, ibd1 = 1, ibd2 = 1;
                A0 = hap_buf[cmp_idx*2]-'0';
                A1 = hap_buf[(cmp_idx*2)+2]-'0';
                
                ibd0 = find_pDgf(f, pDg00, pDg01, pDg11);
                ibd1 = find_pDgIBD1(A0, A1, f, pDg00, pDg01, pDg11);
                if (A0 == 0 && A1 == 0) {
                    ibd2 = pDg00;
                }
                else if ( (A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0) ) {
                    ibd2 = pDg01;
                }
                else if (A0 == 1 && A1 == 1) {
                    ibd2 = pDg11;
                }

                // if no data left after down-sampling,
                // still mark site as processed & include
                // in final coverage distribution/per-site output,
                // but skip when aggregating likelihoods
                if (n_ref+n_alt < 1) {
                    fprintf( tab_fp, "%s\t%s\t%lu\t%s\t%s\t%lf\t%u\t%u\t%u\t%hu\t%hu\t%e\t%e\t%e\n", 
                             pul->chr, id, pos, ref, alt, f, pul->cov, n_ref, n_alt, 
                             A0, A1, ibd0, ibd1, ibd2 ); 
                    processed++;
                    continue;
                }
                
                sum_ibd0 *= ibd0;
                sum_ibd1 *= ibd1; 
                sum_ibd2 *= ibd2;

                if (LD_MODE) {
                    // compare data to each sample in ref panel,
                    // skip if individual matches Pileup name provided via -N
                    // (to avoid including self in background calculation)
                    for (int n = 0; n < data->n_refids; n++) {
                        double ref_ibd2;
                        int ref_idx = data->refids[n].idx;
                        A0 = hap_buf[ref_idx*2]-'0';
                        A1 = hap_buf[(ref_idx*2)+2]-'0';
                        if (A0 == 0 && A1 == 0) {
                            ref_ibd2 = pDg00;
                        }
                        else if ( (A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0) ) {
                            ref_ibd2 = pDg01;
                        }
                        else if (A0 == 1 && A1 == 1) {
                            ref_ibd2 = pDg11;
                        }
                        if (ref_idx != pu_idx) {
                            sum_ibd2_ref[n] *= ref_ibd2;
                        }
                    }
                }
                snp_count++;
                previous_pos = pos;
                if (snp_count == 1) {
                    sgmt_start = pos;
                }
                else if (snp_count == WINDOWSIZE) {
                    sgmt_end = previous_pos;
                }
                fprintf( tab_fp, "%s\t%s\t%lu\t%s\t%s\t%lf\t%u\t%u\t%u\t%c\t%c\t%e\t%e\t%e\n", 
                                   pul->chr, id, pos, ref, alt, f, pul->cov, n_ref, n_alt, 
                                   hap_buf[cmp_idx*2], hap_buf[(cmp_idx*2)+2], ibd0, ibd1, ibd2 ); 
                processed++;
            }
            int n_refpanel = data->n_refids; // number of reference samples
            if (LD_MODE) {
                // take average of IBD2 likelihoods over all ref samples
                double sum_ibd0_LD = 0;
                for (int k = 0; k < data->n_refids; k++) {
                    if (data->refids[k].idx != pu_idx) {
                        sum_ibd0_LD += sum_ibd2_ref[k];
                    }
                    else {
                        n_refpanel--;
                    }
                }
                fprintf( sum_fp, "%d\t%lu\t%lu\t%e\t%e\t%e\t%d\n", sgmt_count, sgmt_start,
                         sgmt_end, sum_ibd0_LD/n_refpanel, sum_ibd1, sum_ibd2, snp_count );
            }
            else {
                fprintf( sum_fp, "%d\t%lu\t%lu\t%e\t%e\t%e\t%d\n", sgmt_count, sgmt_start,
                         sgmt_end, sum_ibd0, sum_ibd1, sum_ibd2, snp_count );
            }
            sgmt_count++;
        }
        fprintf( tab_fp, "# FINAL COVERAGE DISTRIBUTION:\n" );
        fprintf( tab_fp, "# COVERAGE N_SITES\n" );
        for (int cov = 0; cov < USER_MAX_COV+1; cov++) {
            fprintf( tab_fp, "# %d %lu\n", cov, final_dist[cov] );
        }
        fprintf( tab_fp, "# FINAL MEAN DEPTH = %lf\n", (double)final_total_cov/processed );
        fprintf( tab_fp , "## Number of sites processed: %lu\n", processed );
        fprintf( tab_fp, "## Number of sites skipped: %lu\n", skipped );
        fclose(tab_fp);
        fclose(sum_fp);
        rewind_FS(hap_fp);
        rewind_FS(legend_fp);
    }
    return 0;
}   


/* Prints help message */
void print_help(int code) {
    fprintf( stderr, "IBDGem-2.0: Compares low-coverage sequencing data from an unknown sample to known genotypes\n" );
    fprintf( stderr, "            from a reference individual/panel and calculates the likelihood that the samples\n" );
    fprintf( stderr, "            share 0, 1, or 2 IBD chromosomes.\n\n" );
    fprintf( stderr, "Usage: ./ibdgem [--LD] -H [hap-file] -L [legend-file] -I [indv-file] -P [pileup-file] [other options...]\n" );
    fprintf( stderr, "       OR ./ibdgem [--LD] -V [vcf-file] -P [pileup-file] [other options...]\n" );
    fprintf( stderr, "--LD                            Linkage disequilibrium mode ON (default: OFF)\n" );
    fprintf( stderr, "-V, --vcf  FILE                 VCF file (required if using VCF)\n" );
    fprintf( stderr, "-H, --hap  FILE                 HAP file (required if using IMPUTE)\n" );
    fprintf( stderr, "-L, --legend  FILE              LEGEND file (required if using IMPUTE)\n" );
    fprintf( stderr, "-I, --indv  FILE                INDV file (required if using IMPUTE)\n" );
    fprintf( stderr, "-P, --pileup  FILE              PILEUP file (required)\n" );
    fprintf( stderr, "-N, --pileup-name  STR          Name of Pileup sample (default: UNKWN)\n" );
    fprintf( stderr, "-A, --allele-freqs  FILE        File containing allele frequencies from a background panel;\n" );
    fprintf( stderr, "                                   must be sorted & whitespace-delimited with columns CHROM, POS, AF;\n" );
    fprintf( stderr, "                                   use in conjunction with -c if includes multiple chromosomes\n" );
    fprintf( stderr, "                                   (default: calculate AF from input genotypes)\n" );
    fprintf( stderr, "-S, --sample-list  FILE         File containing subset of samples to compare the\n" );
    fprintf( stderr, "                                   sequencing data against; one line per sample\n" );
    fprintf( stderr, "                                   (default: compare against all samples in genotype file)\n" );
    fprintf( stderr, "-s, --sample  STR               Sample(s) to compare the sequencing data against; comma-separated\n" );
    fprintf( stderr, "                                   without spaces if more than one (e.g. sample1,sample2,etc.)\n" );
    fprintf( stderr, "-B, --background-list  FILE     File containing subset of samples to be used as the background panel\n" );
    fprintf( stderr, "                                   for calculating IBD0 in LD mode; one line per sample\n" );
    fprintf( stderr, "                                   (default: use all samples in genotype file as background)\n" );
    fprintf( stderr, "-p, --positions  FILE           List of sites to compare; can be in position list format with 2 columns\n" );
    fprintf( stderr, "                                   CHROM, POS (1-based coordinates) or BED format (0-based coordinates);\n" );
    fprintf( stderr, "                                   use in conjunction with -c if includes multiple chromosomes\n" );
    fprintf( stderr, "                                   (default: perform comparison at all sites)\n" );
    fprintf( stderr, "-q, --min-qual  FLOAT           Genotype quality minimum when using VCF input (default: no minimum)\n" );
    fprintf( stderr, "-M, --max-cov  INT              Maximum estimated coverage of Pileup data (default: 20)\n" );
    fprintf( stderr, "-F, --max-af  FLOAT             Maximum alternate allele frequency (default: 1)\n" );
    fprintf( stderr, "-f, --min-af  FLOAT             Minimum alternate allele frequency (default: 0)\n" );
    fprintf( stderr, "-D, --downsample-cov  FLOAT     Down-sample to this fold-coverage depth\n" );
    fprintf( stderr, "-w, --window-size  INT          Number of sites per genomic segment over which likelihood results\n" ); 
    fprintf( stderr, "                                   are summarized/aggregated (default: 100)\n" );
    fprintf( stderr, "-O, --out-dir  STR              Path to output directory (default: output to current directory)\n" );
    fprintf( stderr, "-c, --chromosome  STR           Chromosome on which the comparison is done; if not specified,\n" );
    fprintf( stderr, "                                   will assume that all inputs are on one single chromosome\n" );
    fprintf( stderr, "-e, --error-rate  FLOAT         Error rate of sequencing platform (default: 0.02)\n" );
    fprintf( stderr, "-v, --variable-sites-only       If set, make output only for sites that are not\n" );
    fprintf( stderr, "                                   homozygous reference in the genotype file for this sample\n" );
    fprintf( stderr, "-h, --help                      Show this help message and exit\n\n" );
    fprintf( stderr, "Format of likelihood table is tab-delimited with columns:\n" );
    fprintf( stderr, "CHR, rsID, POS, REF, ALT, AF, DP, SQ_NREF, SQ_NALT, GT_A0, GT_A1, LIBD0, LIBD1, LIBD2\n\n" );
    fprintf( stderr, "Format of summary file is tab-delimited with columns:\n" );
    fprintf( stderr, "SEGMENT, START, END, LIBD0, LIBD1, LIBD2, NUM_SITES\n" );
    exit(code);
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {

    clock_t start, end;
    double time_elapsed;
    start = clock();

    int option;
    char* opts = ":V:H:L:I:P:w:N:A:S:s:B:p:q:M:F:f:D:O:c:e:vh";

    char vcf_fn[MAX_FN_LEN];
    char hap_fn[MAX_FN_LEN];
    char legend_fn[MAX_FN_LEN];
    char indv_fn[MAX_FN_LEN];
    char pu_fn[MAX_FN_LEN];
    char sample_fn[MAX_FN_LEN];
    char ref_fn[MAX_FN_LEN];
    char af_fn[MAX_FN_LEN];
    char pos_fn[MAX_FN_LEN];
    char user_cmd[MAX_FIELD_WIDTH];
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
            case 0:
                break;
            case 'V':
                IN_VCF = 1;
                strcpy(vcf_fn, optarg);
                break;
            case 'H':
                IN_IMPUTE = 1;
                strcpy(hap_fn, optarg);
                break;
            case 'L':
                IN_IMPUTE = 1;
                strcpy(legend_fn, optarg);
                break;
            case 'I':
                IN_IMPUTE = 1;
                strcpy(indv_fn, optarg);
                break;
            case 'P':
                strcpy(pu_fn, optarg);
                break;
            case 'w':
                WINDOWSIZE = atoi(optarg);
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
            case 'B':
                OPT_B = 1;
                strcpy(ref_fn, optarg);
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
            case 'q':
                USER_MIN_QUAL = atof(optarg);
                break;
            case 'c':
                uchr = optarg;
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
    if ( OPT_D && target_dp <= 0 ) {
        fprintf( stderr, "[::] ERROR: Invalid down-sample coverage (-D) of %.2f (must be > 0).\n", target_dp );
        exit(0);
    }
    if ( USER_MIN_AF < 0 ) {
        fprintf( stderr, "[::] ERROR: Invalid minimum alternate allele frequency (-f) of %.2f (must be >= 0).\n", USER_MIN_AF );
        exit(0);
    }
    if ( USER_MAX_AF > 1 ) {
        fprintf( stderr, "[::] ERROR: Invalid maximum alternate allele frequency (-F) of %.2f (must be <= 1).\n", USER_MAX_AF );
        exit(0);
    }
    if ( USER_MAX_COV < 1 ) {
        fprintf( stderr, "[::] ERROR: Invalid maximum estimated coverage (-M) of %u (must be >= 1).\n", USER_MAX_COV );
        exit(0);
    }
    if ( USER_MIN_QUAL < 0 ) {
        fprintf( stderr, "[::] ERROR: Invalid genotype quality minimum (-q) of %.2f (must be >= 0).\n", USER_MIN_QUAL );
        exit(0);
    }
    if ( WINDOWSIZE < 2 ) {
        fprintf( stderr, "[::] ERROR: Invalid window size (-w) of %d (must be >= 2).\n", WINDOWSIZE );
        exit(0);
    }

    // first copy the first argument into user_cmd to initialize it
    strcpy(user_cmd, argv[0]);
    strcat(user_cmd, " ");
    for (int i = 1; i < argc; i++) {
        strcat(user_cmd, argv[i]);
        strcat(user_cmd, " ");
    }

    Pu_chr* puc = init_Pu_chr(pu_fn, uchr);
    if (!puc) {
        destroy_Pu_chr(puc);
        fprintf( stderr, "[::] ERROR parsing Pileup data; make sure input is valid.\n" );
        exit(1);
    }

    Comp_dt* data = malloc(sizeof(Comp_dt));
    data->uids = NULL;
    data->refids = NULL;
    data->upos = NULL;
    data->uaf = NULL;
    if (OPT_A) {
        int af_status = read_af(data, af_fn, uchr);
        if (af_status == 1) {
            destroy_CD(data);
            destroy_Pu_chr(puc);
            exit(1);
        }
    }

    if (OPT_P) {
        int pos_status = read_pos(data, pos_fn, uchr);
        if (pos_status == 1) {
            destroy_CD(data);
            destroy_Pu_chr(puc);
            exit(1);
        }
    }

    if (!IN_VCF && !IN_IMPUTE) {
        fprintf( stderr, "[::] ERROR: Missing genotype files.\n" );
        destroy_Pu_chr(puc);
        destroy_CD(data);
        exit(1);
    }

    else if (IN_VCF && IN_IMPUTE) {
        fprintf( stderr, "[::] ERROR: 2 types of genotype inputs detected. Please choose either IMPUTE or VCF format.\n" );
        destroy_Pu_chr(puc);
        destroy_CD(data);
        exit(1);
    }

    else if (IN_VCF) {
        File_Src* vcf_fp = init_FS(vcf_fn);
        if (!vcf_fp) {
            destroy_Pu_chr(puc);
            destroy_CD(data);
            destroy_FS(vcf_fp);
            exit(1);
        }
        
        // skip metadata lines
        char line[MAX_LINE_LEN], header[MAX_LINE_LEN];
        get_line_FS(vcf_fp, line);
        while (strstr(line, "##") == line) {
            get_line_FS(vcf_fp, line);
        }

        if (sscanf(line, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t%30720[^\n]", header) == 1) {
            int sample_flag = vcf_parse_samples(header, data);
            if (sample_flag) {
                destroy_Pu_chr(puc);
                destroy_CD(data);
                destroy_FS(vcf_fp);
                exit(1);
            }
        }
        rewind_FS(vcf_fp);
        
        if (OPT_S1) { 
            int read_samples_res = read_sf(data, sample_fn);
            if (read_samples_res) {
                destroy_CD(data);
                destroy_Pu_chr(puc);
                destroy_FS(vcf_fp);
                exit(1);
            }
        }

        else if (OPT_S2) {
            int read_samples_res = read_scmd(data, sample_str);
            if (read_samples_res) {
                destroy_CD(data);
                destroy_Pu_chr(puc);
                destroy_FS(vcf_fp);
                exit(1);
            }
        }

        if (OPT_B) {
            int read_refpanel_res = read_rf(data, ref_fn);
            if (read_refpanel_res) {
                destroy_CD(data);
                destroy_Pu_chr(puc);
                destroy_FS(vcf_fp);
                exit(1);
            }
        }

        unsigned long** nCk = init_nCk(USER_MAX_COV);
        Dpdist* pu_dist = malloc(sizeof(Dpdist));
        memset(pu_dist->dist, 0, sizeof(unsigned long) * (USER_MAX_COV+1));
        compare_vcf(vcf_fp, data, puc, nCk, target_dp, pu_dist, out_dir, user_cmd);
        destroy_FS(vcf_fp);
        destroy_nCk(nCk, USER_MAX_COV);
        free(pu_dist);
    }

    else if (IN_IMPUTE) {
        File_Src* hap_fp = init_FS(hap_fn);
        File_Src* legend_fp = init_FS(legend_fn);
        if (!hap_fp || !legend_fp) {
            destroy_CD(data);
            destroy_Pu_chr(puc);
            fprintf( stderr, "[::] ERROR parsing hap/legend/indv data; make sure inputs are valid.\n" );
            exit(1);
        }

        int res = read_indv(data, indv_fn);
        if (res) {
            destroy_CD(data);
            destroy_Pu_chr(puc);
            destroy_FS(hap_fp);
            destroy_FS(legend_fp);
            exit(1);
        }

        if (OPT_S1) { 
            int read_samples_res = read_sf(data, sample_fn);
            if (read_samples_res) {
                destroy_CD(data);
                destroy_Pu_chr(puc);
                destroy_FS(hap_fp);
                destroy_FS(legend_fp);
                exit(1);
            }
        }

        else if (OPT_S2) {
            int read_samples_res = read_scmd(data, sample_str);
            if (read_samples_res) {
                destroy_CD(data);
                destroy_Pu_chr(puc);
                destroy_FS(hap_fp);
                destroy_FS(legend_fp);
                exit(1);
            }
        }

        if (OPT_B) {
            int read_refpanel_res = read_rf(data, ref_fn);
            if (read_refpanel_res) {
                destroy_CD(data);
                destroy_Pu_chr(puc);
                destroy_FS(hap_fp);
                destroy_FS(legend_fp);
                exit(1);
            }
        }
        
        unsigned long** nCk = init_nCk(USER_MAX_COV);
        Dpdist* pu_dist = malloc(sizeof(Dpdist));
        memset(pu_dist->dist, 0, sizeof(unsigned long) * (USER_MAX_COV+1));
        compare_impute(hap_fp, legend_fp, data, puc, nCk, target_dp, pu_dist, out_dir, user_cmd);
        destroy_FS(hap_fp);
        destroy_FS(legend_fp);
        destroy_nCk(nCk, USER_MAX_COV);
        free(pu_dist);
    }
    destroy_CD(data);
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
