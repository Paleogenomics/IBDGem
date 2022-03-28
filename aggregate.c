/*** Aggregate: Module for aggregating IBDGem's generated likelihood 
 * calculations over a certain number of sites.
***/

#include "file-io.h"

static unsigned int L_DP = 1;
static unsigned int H_DP = 5;
static unsigned int BSIZE = 100;

int passes_filters(int n_obs) {
    if ( (n_obs < L_DP) || (n_obs > H_DP) ) {
        return 0;
    }
    return 1;
}

int main(int argc, char* argv[]) {

    int option;
    char* opts = ":L:H:n:c:";
    char* comp_fn = NULL;
    int n_sites = 0;
    double af, libd0, libd1, libd2;
    unsigned int dp, a0, a1, n_ref, n_alt;
    size_t pos, bend, bstart = 1;

    while ( (option = getopt(argc, argv, opts)) != -1 ) {
        switch (option) {
            case 'L':
                L_DP = atoi(optarg);
                break;
            case 'H':
                H_DP = atoi(optarg);
                break;
            case 'n':
                BSIZE = atoi(optarg);
                break;
            case 'c':
                comp_fn = strdup(optarg);
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
                fprintf(stderr, "Error parsing command-line options.\n");
                exit(0);
        }
    }
    
    for (int i = optind; i < argc; i++) {
        printf("Non-option argument %s\n", argv[i]);
    }

    if (!comp_fn) {
        fprintf( stderr, "AGGREGATE: Summarizes the output of IBDGem\n" );
	    fprintf( stderr, "Writes the aggregated likelihoods across bins of the input file\n" );
	    fprintf( stderr, "These bins can be plotted to see regional trends\n" );
	    fprintf( stderr, "Input data must be sorted by position\n" );
        fprintf( stderr, "-c (required) <input comparison table from IBDGem>\n" );  
	    fprintf( stderr, "-n            <number of sites per bin (default: 100)>\n" );
	    fprintf( stderr, "-L            <low coverage cutoff for comparison data (default: 1)>\n" );
	    fprintf( stderr, "-H            <high coverage cutoff for comparison data (default: 5)>\n" );
        fprintf( stderr, "Format of output table is tab-delimited with columns:\n" );
        fprintf( stderr, "POS_START POS_END AGGR_IBD0 AGGR_IBD1 AGGR_IBD2 NUM_SITES\n" );
        exit(0);
    }

    File_Src* f = init_FS(comp_fn);
    if (!f) {
        exit(1);
    }
    long double aggr_ibd0 = 1.0;
    long double aggr_ibd1 = 1.0;
    long double aggr_ibd2 = 1.0;
    char line[MAX_LINE_LEN+1];
    while ( get_line_FS(f, line) ) {
        if ( sscanf(line, "%zu\t%*s\t%*s\t%*s\t%lf\t%u\t%u\t%u\t%u\t%u\t%le\t%le\t%le", 
            &pos, &af, &dp, &a0, &a1, &n_ref, &n_alt, &libd0, &libd1, &libd2) == 10 ) {
                
            if ( passes_filters(n_ref+n_alt) ) {
                if ( n_sites == BSIZE ) {
                    printf( "%zu\t%zu\t%.7Le\t%.7Le\t%.7Le\t%d\n",
                            bstart, bend, aggr_ibd0, aggr_ibd1, aggr_ibd2, n_sites );
                    n_sites = 0;
                    bstart = pos;
                    aggr_ibd0 = 1.0;
                    aggr_ibd1 = 1.0;
                    aggr_ibd2 = 1.0;
                }
                aggr_ibd0 *= libd0;
                aggr_ibd1 *= libd1;
                aggr_ibd2 *= libd2;
                bend = pos;
                n_sites++;
            }
        }
    }
    if (n_sites > 0) {
        printf( "%zu\t%zu\t%.7Le\t%.7Le\t%.7Le\t%d\n",
                bstart, bend, aggr_ibd0, aggr_ibd1, aggr_ibd2, n_sites );
    }
    destroy_FS(f);
    free(comp_fn);
    return 0;
}