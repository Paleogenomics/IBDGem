/***
 * HiddenGem: Module for inferring IBD0, IBD1, and IBD2 segments across a genomic region
 * between 2 samples using aggregated likelihood calculations from IBDGem.
***/

#include <getopt.h>
#include "file-io.h"

#define MAX_SEGMENTS 12288
#define P1 1
#define P2 2
#define P3 3

static double IBD0_IBD1_penalty = 0.001; // assuming 0.01 crossover/1Mbp & ~100kbp/1 bin, 4 * 0.01 recom/meiosis * 500kbp
static double IBD0_IBD2_penalty = 0.000001;
static double IBD1_IBD2_penalty = 0.001;

/** Long options table **/
static struct option longopts[] = {
    { "summary", required_argument, 0, 's' },
    { "chr",     required_argument, 0, 'c' },
    { "out",     required_argument, 0, 'o' },
    { "p01",     required_argument, 0, P1  },
    { "p02",     required_argument, 0, P2  },
    { "p12",     required_argument, 0, P3  },
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

/*** MAIN PROGRAM ***/
int main(int argc, char* argv[]) {

    int option;
    char* opts = ":s:c:o:";
    char* summary_fn = NULL;
    char* out_fn = NULL;
    char* chr = NULL;

    while ( (option = getopt_long(argc, argv, opts, longopts, NULL)) != -1 ) {
        switch (option) {
            case 's':
                summary_fn = strdup(optarg);
                break;
            case 'c':
                chr = strdup(optarg);
                break;
            case 'o':
                out_fn = strdup(optarg);
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
            case ':':
                fprintf( stderr, "Option missing required argument.\n" );
                exit(0);
            case '?':
                if (isprint(optopt)) {
                    fprintf( stderr, "Invalid option -%c.\n", optopt );
                }
                else {
                    fprintf ( stderr, "Invalid option character '\\x%x'.\n", optopt );
                }
                break;
            default:
                fprintf( stderr, "Error parsing command-line options.\n" );
                exit(0);
        }
    }

    for (int i = optind; i < argc; i++) {
            fprintf( stderr, "Non-option argument %s.\n", argv[i] );
    }

    if (!summary_fn) {
        fprintf( stderr, "HiddenGem: Find most probable path of IBD states across a\n" );
        fprintf( stderr, "genomic region using aggregated likelihoods from IBDGEM.\n" );
        fprintf( stderr, "Usage: ./hiddengem [--p01 IBD0_IBD1_PENALTY --p02 IBD0_IBD2_PENALTY --p12 IBD1_IBD2_PENALTY] --summary (-s) FILE\n" );
        fprintf( stderr, "--summary (-s) <file containing aggregated likelihood results from IBDGem>\n" );
        fprintf( stderr, "--p01 <penalty for changing from state IBD0 to IBD1 and vice versa>\n" );
        fprintf( stderr, "--p02 <penalty for changing from state IBD0 to IBD2 and vice versa>\n" );
        fprintf( stderr, "--p12 <penalty for changing from state IBD1 to IBD2 and vice versa>\n" );
        fprintf( stderr, "Format of output table is tab-delimited with columns:\n" );
        fprintf( stderr, "Bin, IBD0_Score, IBD1_Score, IBD2_Score, Inferred_State\n" );
        exit(0);
    }

    Ibd_genome* ig = init_summary(summary_fn);
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

    printf( "%% IBD0 (n = %.0f): %.2f\n", n_ibd0, (n_ibd0/n)*100 );
    printf( "%% IBD1 (n = %.0f): %.2f\n", n_ibd1, (n_ibd1/n)*100 );
    printf( "%% IBD2 (n = %.0f): %.2f\n", n_ibd2, (n_ibd2/n)*100 );
    
    //fprintf(fpout, "# p01=%lf, p02=%lf, p12=%lf\n",
    //        IBD0_IBD1_penalty, IBD0_IBD2_penalty, IBD1_IBD2_penalty);
    //for (int i = 0; i < n; i++) {
    //    fprintf(fpout, "%d\n", ibd_path[i]);
    //}

    free(summary_fn);
    free(out_fn);
    free(chr);
    free(ig);
    destroy_score_tab(score_tab);
    return 0;
}