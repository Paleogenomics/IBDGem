#include "load_i2.h"
#include "file-io.h"

/* Private function for parsing .hap file */
static void read_hap(File_Src* hf, Impute2* i2) {
    
    /* initialize array of int pointers */
    unsigned short** haps = NULL;
    
    /* initialize number of elements (i.e. haplotypes) in each int array */
    int n_haps;

    /* current array index */
    int curr = 0;
    
    int n = 0;
    char line[MAX_LINE_LEN+1];

    while (get_line_FS(hf, line)) {

        /* extend array as needed */
        if (curr == n) {
            n += STEPSIZE;
            unsigned short** tmp = realloc(haps, n * sizeof(unsigned short*));
            if (!tmp) {
                fprintf(stderr, "Can't realloc.\n");
                exit(1);
            }
            haps = tmp;
        }
        int line_len = strlen(line);

        /* number of haplotypes = half of actual line length */
        n_haps = line_len/2;
        
        unsigned short* allele_arr = malloc(n_haps * sizeof(unsigned short));

        /* skip white space - increase index by 2 */
        for (int i = 0; i < line_len; i += 2) {
            char allele = line[i];
            switch(allele) {
                case '0':
                    allele_arr[i/2] = (unsigned short)0;
                    break;
                case '1':
                    allele_arr[i/2] = (unsigned short)1;
                    break;
                case '2':
                    allele_arr[i/2] = (unsigned short)2;
                    break;
                case '3':
                    allele_arr[i/2] = (unsigned short)3;
                    break;
                case '4':
                    allele_arr[i/2] = (unsigned short)4;
                    break;
                default:
                    fprintf(stderr, "Invalid allele (must be 0-4)\n");
                    exit(1);
            }
        }
        haps[curr] = allele_arr;
        curr++;
    }
    i2->haps = haps;
    i2->num_sites = curr;
    i2->num_haps = n_haps;
}

/* Private function for parsing .legend file */
static void read_legend(File_Src* lf, Impute2* i2) {

    /* obtain number of variants */
    int n = i2->num_sites;

    /* number of fields/columns in .legend file */
    int n_cols = 4;

    /* allocate space for each field/column */
    char** ids = malloc(n * sizeof(char*));
    char** pos = malloc(n * sizeof(char*));
    char** ref_alleles = malloc(n * sizeof(char*));
    char** alt_alleles = malloc(n * sizeof(char*));

    char*** all_cols[4] = {&ids, &pos, &ref_alleles, &alt_alleles};
    
    char line[MAX_LINE_LEN+1];

    /* current line index */
    int i = 0;
    
    /* skip the first (header) line */
    get_line_FS(lf, line);

    while (get_line_FS(lf, line)) {
        if (i == n) {
            break;
        }

        /* copy string into buffer (prevent strtok from modifying line) */
        char buf[MAX_LINE_LEN+1];
        strcpy(buf, line);

        char* token = strtok(buf, " ");
        for (int j = 0; j < n_cols; j++) {
            /* make sure the legend file contains all 4 fields */
            if (!token) {
            fprintf(stderr, "Error parsing .legend, check if input file has at least 4 columns (ID, pos, allele0, allele1)\n");
            exit(1);
            }
            /* trim off new line char in the last column */
            if (j == (n_cols-1)) {
                token[strlen(token)-1] = '\0';
            } 
            int tlen = strlen(token);
            char* str = malloc((tlen+1) * sizeof(char));
            strcpy(str, token);
            (*all_cols[j])[i] = str;
            token = strtok(NULL,  " ");
        }
        i++;
    }
    i2->ids = ids;
    i2->pos = pos;
    i2->ref_alleles = ref_alleles;
    i2->alt_alleles = alt_alleles;
}

/* Private function for parsing .indv file */
static void read_indv(File_Src* inf, Impute2* i2) {
    
    /* obtain number of samples */
    int n = (i2->num_haps)/2;

    /* allocate space for array */
    char** samples = malloc(n * sizeof(char*));

    /* current array index */
    int i = 0;

    char line[MAX_LINE_LEN+1];
    
    while (get_line_FS(inf, line)) {
        if (i == n) {
            break;
        }
        /* trim off new line char */
        line[strlen(line)-1] = '\0';
        int line_len = strlen(line);
        char* str = malloc((line_len+1) * sizeof(char));
        strcpy(str, line);
        samples[i] = str;
        i++;
    }
    i2->samples = samples;
}

/* Private function for freeing memory allocated for array of char pointers */
static int free_cptrarr(char** ptrarr, int n) {
    if (!ptrarr) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        free(ptrarr[i]);
    }
    free(ptrarr);
    return 0;
}

/* Private function for freeing memory allocated for array of ushort pointers */
static int free_iptrarr(unsigned short** ptrarr, int n) {
    if (!ptrarr) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        free(ptrarr[i]);
    }
    free(ptrarr);
    return 0;
}


/** MAIN FUNCTIONS **/
Impute2* init_I2(const char* hap_fn, const char* legend_fn, const char* indv_fn) {

    File_Src* hf = init_FS(hap_fn);
    File_Src* lf = init_FS(legend_fn);
    File_Src* inf = init_FS(indv_fn);

    Impute2* i2 = malloc(sizeof(Impute2));
    read_hap(hf, i2);
    read_legend(lf, i2);
    read_indv(inf, i2);

    destroy_FS(hf);
    destroy_FS(lf);
    destroy_FS(inf);
    return i2;
}

int destroy_I2(Impute2* i2) {
    if (!i2) {
        return 0;
    }
    int nsites = i2->num_sites;
    int nsmpls = (i2->num_haps)/2;
    free_iptrarr(i2->haps, nsites);
    free_cptrarr(i2->samples, nsmpls);
    free_cptrarr(i2->ids, nsites);
    free_cptrarr(i2->pos, nsites);
    free_cptrarr(i2->ref_alleles, nsites);
    free_cptrarr(i2->alt_alleles, nsites);

    free(i2);
    return 0;
}
