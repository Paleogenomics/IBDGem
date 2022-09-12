/***
 * Aggregate: Module for aggregating IBDGem's likelihood results over a fixed number of sites.
***/

#include <getopt.h>
#include "file-io.h"

static unsigned int L_DP = 1;
static unsigned int H_DP = 5;
static unsigned int BSIZE = 100;

/** Long options table **/
static struct option longopts[] = {
    { "comparison-file",      required_argument, 0, 'c' },
    { "num-sites",            required_argument, 0, 'n' },
    { "max-cov",              required_argument, 0, 'H' },
    { "min-cov",              required_argument, 0, 'L' },
    { "help",                 no_argument, 0, 'h' },
    { 0, 0, 0, 0}
};


int passes_filters(int n_obs) {
    if ( (n_obs < L_DP) || (n_obs > H_DP) ) {
        return 0;
    }
    return 1;
}


void print_help(int code) {
    fprintf( stderr, "AGGREGATE: Summarizes IBDGem output by partitioning the genomic region into bins containing\n" );
    fprintf( stderr, "           a fixed number of SNPs and calculates the aggregated likelihoods in each bin.\n" );
    fprintf( stderr, "           These bins can be plotted to see regional trends. Input data must be sorted by position.\n\n" );
    fprintf( stderr, "Usage: ./aggregate -c [comparison-file] [other options...] >[out-file]\n" );
    fprintf( stderr, "-c, --comparison-file  FILE      Comparison table from IBDGem (required)\n" ); 
	fprintf( stderr, "-n, --num-sites  INT             Number of sites per bin (default: 100)\n" );
    fprintf( stderr, "-H, --max-cov  INT               High coverage cutoff for comparison data (default: 5)\n" );
	fprintf( stderr, "-L, --min-cov  INT               Low coverage cutoff for comparison data (default: 1)\n" );
	fprintf( stderr, "-h, --help                       Show this help message and exit\n\n" );
    fprintf( stderr, "Format of output table is tab-delimited with columns:\n" );
    fprintf( stderr, "POS_START POS_END AGGR_LIBD0 AGGR_LIBD1 AGGR_LIBD2 NUM_SITES\n" );
    exit(code);
}


int main(int argc, char* argv[]) {

    int option;
    char* opts = ":c:n:H:L:h";
    char comp_fn[MAX_FN_LEN];
    int n_sites = 0;
    double af, libd0, libd1, libd2;
    unsigned int dp, a0, a1, n_ref, n_alt;
    size_t pos, bend, bstart = 1;

    if (argc == 1) {
        print_help(0);
    }

    while ( (option = getopt_long(argc, argv, opts, longopts, NULL)) != -1 ) {
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
                strcpy(comp_fn, optarg);
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
    
    File_Src* f = init_FS(comp_fn);
    if (!f) {
        fprintf( stderr, "[::] ERROR parsing comparison data; make sure input is valid.\n" );
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
    return 0;
}