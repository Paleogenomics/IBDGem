#include "load-i2.h"
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
                fprintf( stderr, "[::] ERROR in read_hap(): Cannot realloc.\n" );
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
        fprintf( stderr, "[::] ERROR in read_legend(): Cannot parse lines from legend file.\n" );
        free_cparr(ids, n);
        free_cparr(ref, n);
        free_cparr(alt, n);
        free(pos);
        return 1;
    }

    else if (i < n) {
        fprintf( stderr, "[::] ERROR in read_legend(): Fewer sites (%d) than indicated by .hap (%d).\n", i, n );
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

    Sampl* samples = malloc(n * sizeof(Sampl));
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
        fprintf( stderr, "[::] ERROR in read_indv(): Cannot parse lines from indv file.\n" );
        free(samples);
        return 1;
    }
    i2->samples = samples;
    return 0;
}


Sampl* find_sample(Impute2* i2, const char* identifier) {
    for (int i = 0; i < (i2->n_haps)/2; i++) {
        if (strcmp(i2->samples[i].name, identifier) == 0) {
            return &i2->samples[i];
        }
    }
    fprintf( stderr, "Sample %s not found in .indv file\n", identifier );
    return NULL;
}


/* Private function for comparing two Freq objects
   based on their pos values 
   Args: const void* freq1 - pointer to 1st Freq object 
         const void* freq2 - pointer to 2nd Freq object 
   Returns: -1 if freq1->pos < freq2->pos;
             1 if freq1->pos > freq2->pos;
             0 if freq1->pos = freq2->pos;
             2 otherwise */
static int cmp_Freq(const void* freq1, const void* freq2) {
    Freq* f1 = (Freq*) freq1;
    Freq* f2 = (Freq*) freq2;
    size_t pos1, pos2;
    pos1 = f1->pos;
    pos2 = f2->pos;
    if ( pos1 < pos2 ) {
        return -1;
    }
    if ( pos1 > pos2 ) {
        return 1;
    }
    if ( pos1 == pos2 ) {
        return 0;
    }
    return 2;
}


Freq* fetch_freq(Impute2* i2, const size_t pos) {
    Freq key;
    Freq* found_fp;
    key.pos = pos;
    found_fp = bsearch(&key, i2->uaf, i2->n_uaf, sizeof(Freq), cmp_Freq);
    return found_fp;
}


double find_f(Impute2* i2, size_t site_idx) {
    unsigned short* allele_arr = i2->haps[site_idx];
    size_t n_alt = 0; 
    for (int i = 0; i < i2->n_haps; i++) {
        if (*(allele_arr+i) == 1) {
            n_alt++;
        }
    } 
    return (double)n_alt / i2->n_haps;
}


int read_sf(Impute2* i2, const char* sample_fn) {
    File_Src* sf = init_FS(sample_fn);
    if (!sf) {
        return -1;
    }
    Sampl* subset = NULL;
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN+1];

    while (get_line_FS(sf, line)) {
        if (curr == n) {
            n += 100;
            Sampl* tmp = realloc(subset, n * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_sf(): Cannot realloc.\n" );
                if (subset) {
                    free(subset);
                }
                return -1;
            }
            subset = tmp;
        }
        line[strlen(line)-1] = '\0';
        Sampl* sp = find_sample(i2, line);
        if (!sp) {
            fprintf( stderr, "Sample %s not found in input panel.\n", line );
            continue;
        }
        subset[curr] = *sp;
        curr++;
    }  
    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_samples(): No matching samples found in %s.\n", sample_fn );
        destroy_FS(sf);
        free(subset);
        return -1;
    }
    free(i2->samples);
    i2->samples = subset;
    destroy_FS(sf);
    return curr;
}


int read_scmd(Impute2* i2, const char* sample_str) {
    Sampl* subset = NULL;
    int curr = 0;
    int n = 0;
    char buf[MAX_LINE_LEN+1];
    strcpy(buf, sample_str);

    char* sample = strtok(buf, ",");
    while(sample) {
        if (curr == n) {
            n += 100;
            Sampl* tmp = realloc(subset, n * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_scmd(): Cannot realloc.\n" );
                if (subset) {
                    free(subset);
                }
                return -1;
            }
            subset = tmp;
        }
        Sampl* sp = find_sample(i2, sample);
        if (!sp) {
            fprintf( stderr, "Sample %s not found in input panel.\n", sample );
            sample = strtok(NULL, ",");
            continue;
        }
        subset[curr] = *sp;
        curr++;
        sample = strtok(NULL, ",");
    }

    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_scmd(): No matching samples found.\n" );
        free(subset);
        return -1;
    }
    free(i2->samples);
    i2->samples = subset;
    return curr;
}


int read_af(Impute2* i2, const char* af_fn, const char* chr) {
    File_Src* ff = init_FS(af_fn);
    if (!ff) {
        return 1;
    }
    Freq* af = NULL;
    char buf_chr[128];
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN];

    while (get_line_FS(ff, line)) {
        if (curr == n) {
            n += 200;
            Freq* tmp = realloc(af, n * sizeof(Freq));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_af(): Cannot realloc.\n" );
                if (af) {
                    free(af);
                }
                return 1;
            }
            af = tmp;
        }
        if (sscanf(line, "%128s %lu %lf", buf_chr, &af[curr].pos, &af[curr].f) == 3) {
            // if site not from specified chromosome, skip it
            if (chr) {
                if (strcmp(buf_chr, chr) == 0) {
                    curr++;
                }
            }
            else {
                curr++;
            }
        }
    }
    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_af(): Cannot parse lines from %s.\n", af_fn );
        free(af);
        return 1;
    }

    i2->uaf = af;
    i2->n_uaf = curr;
    destroy_FS(ff);
    return 0;
}


int read_pos(Impute2* i2, const char* pos_fn, const char* chr) {
    File_Src* pf = init_FS(pos_fn);
    if (!pf) {
        return 1;
    }

    unsigned long* upos = NULL;
    char buf_chr[128];
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN];

    while (get_line_FS(pf, line)) {
        if (curr == n) {
            n += 2000;
            unsigned long* tmp = realloc(upos, n * sizeof(unsigned long));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_pos(): Cannot realloc.\n" );
                if (upos) {
                    free(upos);
                }
                return 1;
            }
            upos = tmp;
        }
        // checks if either in BED or CHROM,POS format
        if ( (sscanf(line, "%128s %*d %lu", buf_chr, &upos[curr]) == 2) ||
             (sscanf(line, "%128s %lu", buf_chr, &upos[curr]) == 2) ) {
            // if site not from specified chromosome, skip it
            if (chr) {
                if (strcmp(buf_chr, chr) == 0) {
                    curr++;
                }
            }
            else {
                curr++;
            }
        }
    }
    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_pos(): Cannot parse lines from %s.\n", pos_fn );
        free(upos);
        return 1;
    }

    i2->upos = upos;
    i2->n_upos = curr;
    destroy_FS(pf);
    return 0;
}


int site_in_upos(Impute2* i2, const size_t pos) {
    for (int i = 0; i < i2->n_upos; i++) {
        if (pos == i2->upos[i]) {
            return 1;
        }
    }
    return 0;
}


/** MAIN FUNCTIONS **/
Impute2* init_I2(const char* hap_fn, const char* legend_fn, const char* indv_fn) {

    File_Src* hf = init_FS(hap_fn);
    File_Src* lf = init_FS(legend_fn);
    File_Src* inf = init_FS(indv_fn);
    if ( !hf || !lf || !inf ) {
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
    i2->uaf = NULL;
    i2->upos = NULL;
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
    free(i2->upos);
    free(i2->samples);
    free(i2->uaf);
    free(i2);
    return 0;
}
