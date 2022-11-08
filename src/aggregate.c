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