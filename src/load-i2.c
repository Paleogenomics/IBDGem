#include "load-i2.h"
#include "file-io.h"

static int VCF_HDR_PRESENT = 0;
static int VCF_FRMT_PRESENT = 0;

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


/* Private function for determining if site is biallelic 
   Args: const char* alt - pointer to alternate allele string 
   Returns: 1 if biallelic, 0 otherwise */
static int is_biallelic(const char* alt) {
    // see if alt contains a comma
    if (strchr(alt, ',') == NULL) {
        return 1;
    }
    return 0;
}


/* Private function for determining if genotype field is valid
   (first 3 chars are 0s and 1s separated by | or /) 
   Args: const char* gt - pointer to genotype string
         regex_t regex - regular expression for pattern checking 
   Returns: 1 if genotype valid, 0 otherwise */
static int genotype_ok(const char* gt, regex_t regex) {
    int rc = regexec(&regex, gt, 0, NULL, 0);
    if (!rc) {
        return 1;
    }
    return 0;
}


/* Private function for extracting alleles from multiple
   genotype fields
   Args: char* gt - pointer to genotype fields as a single string
         regex_t regex - regular expression for pattern checking
         int n_samples - number of samples
         unsigned short* alleles - pointer to alleles array
   Returns: 0 if successfully parsed, 1 otherwise */
static int vcf_parse_gt(char* gt, regex_t regex, int n_samples, unsigned short* alleles) {
    char* token;
    token = strtok(gt, "\t");
    // skip FORMAT column if present
    if (VCF_FRMT_PRESENT == 1) {
        token = strtok(NULL, "\t");
    }
    for (int i = 0; i < n_samples; i++) {
        if ( genotype_ok(token, regex) )  {
            alleles[i*2] = token[0] - '0';
            alleles[i*2+1] = token[2] - '0';
        }
        else {
            return 1;
        }
        token = strtok(NULL, "\t");
    }

    return 0;
}


/* Private function for parsing sample names from VCF header
   Args: const char* header - pointer to header string
         Impute2* i2 - pointer to Impute genotype info
   Returns: 0 if successfully parsed, 1 otherwise */
static int vcf_parse_samples(const char* header, Impute2* i2) {
    int n = 0;
    int arr_size = 200;
    Sampl* samples = malloc(arr_size * sizeof(Sampl));
    char* token;
    char buf[MAX_LINE_LEN];
    strcpy(buf, header);
    
    token = strtok(buf, "\t");
    // skip FORMAT column if present
    if (strcmp(token, "FORMAT") == 0) {
        VCF_FRMT_PRESENT = 1;
        token = strtok(NULL, "\t");
    }
    while (token != NULL) {
        if (n == arr_size) {
            arr_size += 200;
            Sampl* tmp = realloc(samples, arr_size * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in vcf_parse_samples(): Cannot realloc.\n" );
                free(samples);
                return 1;
            }
            samples = tmp;
        }
        strcpy(samples[n].name, token);
        samples[n].idx = n*2;
        token = strtok(NULL, "\t");
        n++;
    }

    if (n == 0) {
        fprintf( stderr, "[::] ERROR: No samples found.\n" );
        free(samples);
        return 1;
    }
    i2->samples = samples;
    i2->n_haps = n*2;
    return 0;
}


/* Private function for parsing VCF file
   Args: File_Src* vcf - pointer to VCF file
         Impute2* i2 - pointer to Impute genotype info
   Returns: 0 if file read successfully, 1 otherwise */
static int read_vcf(File_Src* vcf, Impute2* i2) {
    char line[MAX_LINE_LEN];
    char header[MAX_LINE_LEN];
    char buf_gt[MAX_LINE_LEN];
    char buf_id[512], buf_ref[512], buf_alt[512], buf_qual[512], buf_fltr[512], buf_info[512];
    unsigned long buf_pos;

    char** ids = NULL, **ref = NULL, **alt = NULL;
    unsigned long* pos = NULL;
    double* qual = NULL;
    unsigned short** haps = NULL;
    int n = 0;
    int arr_size = 0;
    
    // skip metadata lines
    get_line_FS(vcf, line);
    while (strstr(line, "##") == line) {
        get_line_FS(vcf, line);
    }

    if (sscanf(line, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t%30720[^\n]", header) == 1) {
        VCF_HDR_PRESENT = 1;
        int sample_flag = vcf_parse_samples(header, i2);
        if (sample_flag) {
            return 1;
        }
    }
    
    if (!VCF_HDR_PRESENT) {
        fprintf( stderr, "[::] ERROR in read_vcf(): No header found.\n" );
        return 1;           
    }
    
    regex_t regex;
    int rc_flag = regcomp(&regex, "^[01][/|][01].*$", 0);
    if (rc_flag) {
        fprintf( stderr, "[::] ERROR in read_vcf(): Cannot compile regex to check genotype.\n" );
        return 1;
    }

    while (get_line_FS(vcf, line)) {
        if (sscanf(line, "%*s %lu %512[^\t] %512[^\t] %512[^\t] %512[^\t] %512[^\t] %512[^\t] %30720[^\n]",
                   &buf_pos, buf_id, buf_ref, buf_alt, buf_qual, buf_fltr, buf_info, buf_gt) == 8) {

            // skip site if not biallelic        
            if (!is_biallelic(buf_alt)) {
                continue;
            }
            if (n == arr_size) {
                arr_size += STEPSIZE;
                char** tmp_ids = realloc(ids, arr_size * sizeof(char*));
                char** tmp_ref = realloc(ref, arr_size * sizeof(char*));
                char** tmp_alt = realloc(alt, arr_size * sizeof(char*));
                unsigned long* tmp_pos = realloc(pos, arr_size * sizeof(unsigned long));
                double* tmp_qual = realloc(qual, arr_size * sizeof(double));
                unsigned short** tmp_haps = realloc(haps, arr_size * sizeof(unsigned short*));
                if (!tmp_ids || !tmp_ref || !tmp_alt || !tmp_pos || !tmp_haps) {
                    fprintf( stderr, "[::] ERROR in read_vcf(): Cannot realloc.\n" );
		            free_cparr(ids, n);
                    free_cparr(ref, n);
                    free_cparr(alt, n);
                    free_iparr(haps, n);
                    free(pos);
                    free(qual);
                    regfree(&regex);
                    return 1;
                }
                ids = tmp_ids;
                ref = tmp_ref;
                alt = tmp_alt;
                pos = tmp_pos;
                qual = tmp_qual;
                haps = tmp_haps;
            }
            ids[n] = strdup(buf_id);
            ref[n] = strdup(buf_ref);
            alt[n] = strdup(buf_alt);
            pos[n] = buf_pos;
            qual[n] = atof(buf_qual);

            unsigned short* alleles = malloc(i2->n_haps * sizeof(unsigned short));
            int gt_flag = vcf_parse_gt(buf_gt, regex, i2->n_haps/2, alleles);
            if (gt_flag) {
                fprintf( stderr, "Failed to parse genotype fields at %lu. Skipping to next site.\n", pos[n]);
                continue;
            }
            haps[n] = alleles;
            n++;
        }
    }
    
    if (n == 0) {
        fprintf( stderr, "[::] ERROR in read_vcf(): No variants parsed from VCF.\n" );
        regfree(&regex);
        return 1;
    }
    
    i2->n_sites = n;
    i2->ids = ids;
    i2->ref = ref;
    i2->alt = alt;
    i2->pos = pos;
    i2->qual = qual;
    i2->haps = haps;
    regfree(&regex);
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
    int n = 0;
    
    int arr_size = 0;
    char line[MAX_LINE_LEN+1];

    while (get_line_FS(hf, line)) {

        // extend array as needed
        if (n == arr_size) {
            arr_size += STEPSIZE;
            unsigned short** tmp = realloc(haps, arr_size * sizeof(unsigned short*));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_hap(): Cannot realloc.\n" );
                free_iparr(haps, n);
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
                default:
                    fprintf( stderr, "Invalid allele (expect 0 and 1).\n" );
            }
        }
        haps[n] = alleles;
        n++;
    }

    i2->haps = haps;
    i2->n_sites = n;
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
        fprintf( stderr, "[::] ERROR in read_legend(): Failed to parse legend file.\n" );
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
        fprintf( stderr, "[::] ERROR in read_indv(): No samples found.\n" );
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
        fprintf( stderr, "[::] ERROR in read_af(): Failed to parse %s.\n", af_fn );
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
        fprintf( stderr, "[::] ERROR in read_pos(): Failed to parse %s.\n", pos_fn );
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


Impute2* init_vcf(const char* vcf_fn) {
    File_Src* vcf = init_FS(vcf_fn);
    if (!vcf) {
        return NULL;
    }
    Impute2* i2 = malloc(sizeof(Impute2));
    int vcfc = read_vcf(vcf, i2);
    if (vcfc == 1) {
        destroy_FS(vcf);
        free(i2);
        return NULL;
    }
    i2->uaf = NULL;
    i2->upos = NULL;
    destroy_FS(vcf);
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
    free(i2->qual);
    free(i2->upos);
    free(i2->samples);
    free(i2->uaf);
    free(i2);
    return 0;
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