#include "load_i2.h"
#include "file-io.h"


/* Private function for freeing memory allocated by array of char pointers 
   Args: char** ptrarr - pointer to char pointer array 
         int n - number of elements in pointer array 
   Returns: 0 if freed successfully */
static int free_cparr(char** ptrarr, int n) {
    if (!ptrarr) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        free(ptrarr[i]);
    }
    free(ptrarr);
    return 0;
}


/* Private function for freeing memory allocated by array of ushort pointers 
   Args: unsigned short** ptrarr - pointer to ushort pointer array 
         int n - number of elements in pointer array 
   Returns: 0 if freed successfully */
static int free_iparr(unsigned short** ptrarr, int n) {
    if (!ptrarr) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        free(ptrarr[i]);
    }
    free(ptrarr);
    return 0;
}


/* Private function for parsing .hap file 
   Args: File_Src* hf - pointer to .hap file
         Impute2* i2 - pointer to Impute genotype info
   Returns: 0 if file read successfully, 1 otherwise */
static int read_hap(File_Src* hf, Impute2* i2) {
    
    // declare array of int pointers
    unsigned short** haps = NULL;
    
    // declare number of elements (i.e. haplotypes) in each int array
    int n_haps;
    
    // current array index
    int curr = 0;
    
    int n = 0;
    char line[MAX_LINE_LEN+1];

    while (get_line_FS(hf, line)) {

        // extend array as needed
        if (curr == n) {
            n += STEPSIZE;
            unsigned short** tmp = realloc(haps, n * sizeof(unsigned short*));
            if (!tmp) {
                fprintf( stderr, "Error from read_hap(): Cannot realloc.\n" );
                if (haps) {
                    free_iparr(haps, curr);
                }
                return 1;
            }
            haps = tmp;
        }
        int len = strlen(line);

        // number of haplotypes = half of actual line length
        n_haps = len/2;
        
        unsigned short* alleles = malloc(n_haps * sizeof(unsigned short));

        // skip white space -> increase index by 2
        for (int i = 0; i < len; i += 2) {
            char allele = line[i];
            switch(allele) {
                case '0':
                    alleles[i/2] = (unsigned short)0;
                    break;
                case '1':
                    alleles[i/2] = (unsigned short)1;
                    break;
                case '2':
                    alleles[i/2] = (unsigned short)2;
                    break;
                case '3':
                    alleles[i/2] = (unsigned short)3;
                    break;
                case '4':
                    alleles[i/2] = (unsigned short)4;
                    break;
                default:
                    fprintf( stderr, "Invalid allele (expect 0-4).\n" );
            }
        }
        haps[curr] = alleles;
        curr++;
    }

    i2->haps = haps;
    i2->n_sites = curr;
    i2->n_haps = n_haps;
    return 0;
}


/* Private function for parsing .legend file 
   Args: File_Src* lf - pointer to .legend file
         Impute2* i2 - pointer to Impute genotype info
   Returns: 0 if file read successfully, 1 otherwise */
static int read_legend(File_Src* lf, Impute2* i2) {

    int n = i2->n_sites;

    char** ids = malloc(n * sizeof(char*));
    char** ref = malloc(n * sizeof(char*));
    char** alt = malloc(n * sizeof(char*));
    unsigned long* pos = malloc(n * sizeof(unsigned long));
    
    char line[MAX_LINE_LEN+1];
    char* tmp_id = NULL, *tmp_ref = NULL, *tmp_alt = NULL;
    
    int i = 0;
    while (get_line_FS(lf, line)) {
        if (i == n) {
            break;
        }
        if (sscanf(line, "%ms %lu %ms %ms",
            &tmp_id, &pos[i], &tmp_ref, &tmp_alt) == 4) {

            ids[i] = strdup(tmp_id);
            ref[i] = strdup(tmp_ref);
            alt[i] = strdup(tmp_alt);
            i++;
        }
        free(tmp_id);
        free(tmp_ref);
        free(tmp_alt);
    }
    
    if (i == 0) {
        fprintf( stderr, "Error from read_legend(): No valid sites found; make sure .legend file\n" );
        fprintf( stderr, "has at least 4 columns (ID, pos, allele0, allele1).\n" );
        free_cparr(ids, n);
        free_cparr(ref, n);
        free_cparr(alt, n);
        free(pos);
        return 1;
    }

    else if (i < n) {
        fprintf( stderr, "Error from read_legend(): Legend has fewer sites (%d) than indicated by .hap (%d).\n", i, n );
        free_cparr(ids, n);
        free_cparr(ref, n);
        free_cparr(alt, n);
        free(pos);
        return 1;
    }

    i2->ids = ids;
    i2->pos = pos;
    i2->ref = ref;
    i2->alt = alt;
    return 0;
}


/* Private function for parsing .indv file
   Args: File_Src* inf - pointer to .indv file
         Impute2* i2 - pointer to Impute genotype info
   Returns: 0 if file read successfully, 1 otherwise */
static int read_indv(File_Src* inf, Impute2* i2) {
    
    // get number of samples
    int n = (i2->n_haps)/2;

    Sample_info* samples = malloc(n * sizeof(Sample_info));
    char line[MAX_LINE_LEN+1];
    
    int i = 0;
    while (get_line_FS(inf, line)) {
        if (i == n) {
            break;
        }
        
        // trim off new line char
        line[strlen(line)-1] = '\0';
        strcpy(samples[i].name, line);
        samples[i].idx = i*2;
        i++;
    }
    if (i == 0) {
        fprintf( stderr, "Error from read_indv(): No valid samples found.\n" );
        free(samples);
        return 1;
    }
    i2->samples = samples;
    return 0;
}


/** MAIN FUNCTIONS **/
Impute2* init_I2(const char* hap_fn, const char* legend_fn, const char* indv_fn) {

    File_Src* hf = init_FS(hap_fn);
    File_Src* lf = init_FS(legend_fn);
    File_Src* inf = init_FS(indv_fn);
    if ( !hf|| !lf || !inf ) {
        return NULL;
    }
    Impute2* i2 = malloc(sizeof(Impute2));
    int hc = read_hap(hf, i2);
    int lc = read_legend(lf, i2);
    int ic = read_indv(inf, i2);
    if ( hc == 1 || lc == 1 || ic == 1 ) {
        destroy_FS(hf);
        destroy_FS(lf);
        destroy_FS(inf);
        free(i2);
        return NULL;
    }
    destroy_FS(hf);
    destroy_FS(lf);
    destroy_FS(inf);
    return i2;
}


int destroy_I2(Impute2* i2) {
    if (!i2) {
        return 0;
    }
    int n = i2->n_sites;
    free_iparr(i2->haps, n);
    free_cparr(i2->ids, n);
    free_cparr(i2->ref, n);
    free_cparr(i2->alt, n);
    free(i2->pos);
    free(i2->samples);
    free(i2);
    return 0;
}
