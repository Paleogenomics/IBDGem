/***
 * HiddenGem: Module for inferring IBD0, IBD1, and IBD2 segments across a genomic region.
***/

#include <getopt.h>
#include "file-io.h"

#define MAX_SEGMENTS 12288
#define P1 1
#define P2 2
#define P3 3

static double IBD0_IBD1_penalty = 0.001; // assuming 0.01 crossover/1Mbp & ~100kbp/1 bin
static double IBD0_IBD2_penalty = 0.000001;
static double IBD1_IBD2_penalty = 0.001;

/** Long options table **/
static struct option longopts[] = {
    { "summary", required_argument, 0, 's' },
    { "p01",     required_argument, 0, P1  },
    { "p02",     required_argument, 0, P2  }, 
    { "p12",     required_argument, 0, P3  },
    { "help",    no_argument, 0, 'h' },
    { 0, 0, 0, 0}
};

/** Structure for storing data from summary file **/
typedef struct ibd_genome {
    size_t start[MAX_SEGMENTS];
    size_t end[MAX_SEGMENTS];
    double nrm_lhs[3][MAX_SEGMENTS]; // *normalized* likelihoods from IBDGem
    unsigned int n_sites[MAX_SEGMENTS];
    int n_bins;
} Ibd_genome;


/** Structure containing likelihood score info at a
    genomic segment/bin for a specific IBD state **/
typedef struct score {
    int curr_state; // current IBD state
    // maximum total likelihood at this bin:
    // max( likelihood[IBD0/IBD1/IBD2][i-1] + likelihood[curr_state][i] - penalty_if_applicable )
    long double tot_score;
    struct score* previous; // pointer to previous score that gives max likelihood
} Score;


/** Stores data from summary file
 * Args: const char* summary_fn - summary filename
 * Returns: Ibd_genome object with summary data **/
Ibd_genome* init_summary(const char* summary_fn) {
    File_Src* sf = init_FS(summary_fn);
    if (!sf) {
        exit(1);
    }
    Ibd_genome* ig = malloc(sizeof(Ibd_genome));
    double l0, l1, l2;
    char line[MAX_LINE_LEN + 1];
    int i = 0;
    while (get_line_FS(sf, line)) {
        if ( sscanf(line, "%zu\t%zu\t%lf\t%lf\t%lf\t%d",
                &ig->start[i],
                &ig->end[i],
                &l0, &l1, &l2,
                &ig->n_sites[i]) == 6 ) {
            
            // normalize the aggregated likelihood (i.e. convert to probability)
            // for each bin
            ig->nrm_lhs[0][i] = l0/(l0+l1+l2);
            ig->nrm_lhs[1][i] = l1/(l0+l1+l2);
            ig->nrm_lhs[2][i] = l2/(l0+l1+l2);
            i++;
        }
    }
    ig->n_bins = i;
    destroy_FS(sf);
    return ig;
}


/** Finds the index of the largest number in an array
 * Args: long double* arr - pointer to number array
 *       size_t n -  number of elements in array 
 * Returns: index of largest number **/
int find_max_idx(long double* arr, size_t n) {
    int max_idx = 0;
    for (int i = 0; i < n; i++) {
        if (arr[i] > arr[max_idx]) {
            max_idx = i;
        }
    }
    return max_idx;
}

 
/** Populates the score matrix with maximum likelihood scores for every
 * state at each segment/bin; each cell is back-pointed to the previous 
 * state that gives it the highest score
 * Args: Ibd_genome* ig - pointer to summary data
 *       Score** score_tab - pointer to score matrix
 * Returns: 0 if all scores added successfully **/
int calc_score(Ibd_genome* ig, Score** score_tab) {
    int n = ig->n_bins;
    long double from_ibd0, from_ibd1, from_ibd2;
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            for (int s = 0; s < 3; s++) {
                score_tab[s][0].curr_state = s;
                score_tab[s][0].tot_score = ig->nrm_lhs[s][0];
                score_tab[s][0].previous = &score_tab[s][0];
            }
        }
        
        else {
            int max_state;
            for (int s = 0; s < 3; s++) {
                if (s == 0) {
                    from_ibd0 = score_tab[0][i-1].tot_score * ig->nrm_lhs[0][i];
                    from_ibd1 = score_tab[1][i-1].tot_score * ig->nrm_lhs[0][i] * IBD0_IBD1_penalty;
                    from_ibd2 = score_tab[2][i-1].tot_score * ig->nrm_lhs[0][i] * IBD0_IBD2_penalty;
                }
                else if (s == 1) {
                    from_ibd0 = score_tab[0][i-1].tot_score * ig->nrm_lhs[1][i] * IBD0_IBD1_penalty;
                    from_ibd1 = score_tab[1][i-1].tot_score * ig->nrm_lhs[1][i];
                    from_ibd2 = score_tab[2][i-1].tot_score * ig->nrm_lhs[1][i] * IBD1_IBD2_penalty;
                }
                else if (s == 2) {
                    from_ibd0 = score_tab[0][i-1].tot_score * ig->nrm_lhs[2][i] * IBD0_IBD2_penalty;
                    from_ibd1 = score_tab[1][i-1].tot_score * ig->nrm_lhs[2][i] * IBD1_IBD2_penalty;
                    from_ibd2 = score_tab[2][i-1].tot_score * ig->nrm_lhs[2][i];
                }
                long double scores[3] = {from_ibd0, from_ibd1, from_ibd2};
                max_state = find_max_idx(scores, 3);
                score_tab[s][i].curr_state = s;
                score_tab[s][i].tot_score = scores[max_state];
                score_tab[s][i].previous = &score_tab[max_state][i-1];
            }
        }
    }
    return 0;
}


/** Initializes the score matrix
 * Args: int n - number of genomic segments
 * Returns: the score matrix **/
Score** init_score_tab(const int n) {
    Score** score_tab = malloc(3 * sizeof(Score*)); // 3 total IBD states
    for (int i = 0; i < 3; i++) {
        score_tab[i] = malloc(n * sizeof(Score));
    }
    return score_tab;
}


/** Free memory allocated by init_score_tab() **/
int destroy_score_tab(Score** score_tab) {
    if (!score_tab) {
        return 0;
    }
    for (int i = 0; i < 3; i++) {
        free(score_tab[i]);
    }
    free(score_tab);
    return 0;
}


void print_help(int code) {
    fprintf( stderr, "HIDDENGEM: Finds most probable path of IBD states across genomic bins from IBDGem-derived\n" );
    fprintf( stderr, "           aggregated likelihoods.\n\n" );
    fprintf( stderr, "Usage: ./hiddengem -s [summary-file] [other options...] >[out-file]\n" );
    fprintf( stderr, "--summary, -s  FILE      File containing aggregated likelihood results from IBDGem (required)\n" );
    fprintf( stderr, "--p01  FLOAT             Penalty for switching between states IBD0 and IBD1 (default: 1e-3)\n" );
    fprintf( stderr, "--p02  FLOAT             Penalty for switching between states IBD0 and IBD2 (default: 1e-6)\n" );
    fprintf( stderr, "--p12  FLOAT             Penalty for switching between states IBD1 and IBD2 (default: 1e-3)\n" );
    fprintf( stderr, "--help                   Show this help message and exit\n\n" );
    fprintf( stderr, "Format of output table is tab-delimited with columns:\n" );
    fprintf( stderr, "Bin, IBD0_Score, IBD1_Score, IBD2_Score, Inferred_State\n" );
    exit(code);
}


/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {

    int option;
    char* opts = ":s:h";
    char summary_fn[MAX_FN_LEN];

    if (argc == 1) {
        print_help(0);
    }

    while ( (option = getopt_long(argc, argv, opts, longopts, NULL)) != -1 ) {
        switch (option) {
            case 's':
                strcpy(summary_fn, optarg);
                break;
            case P1:
                IBD0_IBD1_penalty = atof(optarg);
                break;
            case P2:
                IBD0_IBD2_penalty = atof(optarg);
                break;
            case P3:
                IBD1_IBD2_penalty = atof(optarg);
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
                    fprintf ( stderr, "Invalid option character.\n" );
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

    Ibd_genome* ig = init_summary(summary_fn);
    if (!ig) {
        fprintf( stderr, "[::] ERROR parsing likelihood data; make sure input is valid.\n" );
        exit(1);
    }
    int n = ig->n_bins;
    Score** score_tab = init_score_tab(n);    
    calc_score(ig, score_tab);
      
    long double final_scores[3] = {score_tab[0][n-1].tot_score,
                                   score_tab[1][n-1].tot_score,
                                   score_tab[2][n-1].tot_score};
    int final_state = find_max_idx(final_scores, 3);
    
    // final state = state that gives maximum score in final bin
    int ibd_path[n];
    ibd_path[n-1] = final_state;
    for (int i = n-2; i >= 0; i--) {
        int s = ibd_path[i+1];
        ibd_path[i] = score_tab[s][i+1].previous->curr_state;
    }

    printf("Bin\tIBD0_Score\tIBD1_Score\tIBD2_Score\tInferred_State\n");
    // bin counts for each state
    double n_ibd0 = 0;
    double n_ibd1 = 0;
    double n_ibd2 = 0;
    for (int i = 0; i < n; i++) {
        if (ibd_path[i] == 0) {
            n_ibd0++;
        }
        else if (ibd_path[i] == 1) {
            n_ibd1++;
        }
        else if (ibd_path[i] == 2) {
            n_ibd2++;
        }
        
        printf("%d\t%.5Le\t%.5Le\t%.5Le\t%d\n", i+1, 
               score_tab[0][i].tot_score,
               score_tab[1][i].tot_score,
               score_tab[2][i].tot_score,
               ibd_path[i]);
    }
    printf( "#%% IBD0 (n = %.0f): %.2f\n", n_ibd0, (n_ibd0/n)*100 );
    printf( "#%% IBD1 (n = %.0f): %.2f\n", n_ibd1, (n_ibd1/n)*100 );
    printf( "#%% IBD2 (n = %.0f): %.2f\n", n_ibd2, (n_ibd2/n)*100 );

    free(ig);
    destroy_score_tab(score_tab);
    return 0;
}



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