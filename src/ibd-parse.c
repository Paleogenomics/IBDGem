#include "ibd-parse.h"
#include "file-io.h"

int read_indv(Comp_dt* data, const char* indv_fn) {
    File_Src* indv_fp = init_FS(indv_fn);
    if (!indv_fp) {
        return 1;
    }
    int n = 0;
    int arr_size = 200;
    Sampl* samples = malloc(arr_size * sizeof(Sampl));
    char line[MAX_LINE_LEN];
    
    while (get_line_FS(indv_fp, line)) {
        if (n == arr_size) {
            arr_size += 200;
            Sampl* tmp = realloc(samples, arr_size * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in parse_indv(): Cannot realloc.\n" );
                free(samples);
                destroy_FS(indv_fp);
                return 1;
            }
            samples = tmp;
        }
        // trim off new line char
        line[strlen(line)-1] = '\0';
        strcpy(samples[n].name, line);
        samples[n].idx = n*2;
        n++;
    }
    if (n == 0) {
        fprintf( stderr, "[::] ERROR: No samples found in .indv file.\n" );
        free(samples);
        destroy_FS(indv_fp);
        return 1;
    }
    data->ids = samples;
    data->n_ids = n;
    destroy_FS(indv_fp);
    return 0;
}


Sampl* find_sample(Comp_dt* data, const char* identifier) {
    for (int i = 0; i < data->n_ids; i++) {
        if (strcmp(data->ids[i].name, identifier) == 0) {
            return &data->ids[i];
        }
    }
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


Freq* fetch_freq(Comp_dt* data, const size_t pos) {
    Freq key;
    Freq* found_fp;
    key.pos = pos;
    found_fp = bsearch(&key, data->uaf, data->n_uaf, sizeof(Freq), cmp_Freq);
    return found_fp;
}


double find_f_impute(const char* hap_str, int n_samples) {
    double n_alt = 0; 
    for (int i = 0; i < (n_samples*4)-1; i+=2) {
        if (hap_str[i] == '1') {
            n_alt++;
        }
    } 
    return n_alt / (n_samples*2);
}


double find_f_vcf(unsigned short* alleles, int n_samples) {
    double n_alt = 0;
    for (int i = 0; i < n_samples*2; i++) {
        if (alleles[i] == 1) {
            n_alt++;
        }
    }
    return n_alt / (n_samples*2);
}


int vcf_parse_samples(const char* header, Comp_dt* data) {
    int n = 0;
    int arr_size = 200;
    Sampl* samples = malloc(arr_size * sizeof(Sampl));
    char* token;
    char buf[MAX_LINE_LEN];
    strcpy(buf, header);
    
    token = strtok(buf, "\t");
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
    data->ids = samples;
    data->n_ids = n;
    return 0;
}


int genotype_ok(const char* gt, regex_t regex) {
    int rc = regexec(&regex, gt, 0, NULL, 0);
    if (!rc) {
        return 1;
    }
    return 0;
}


int vcf_parse_gt(char* gt, regex_t regex, int n_samples, unsigned short* alleles) {
    char* token;
    token = strtok(gt, "\t");
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


int read_sf(Comp_dt* data, const char* sample_fn) {
    File_Src* sample_fp = init_FS(sample_fn);
    if (!sample_fp) {
        return 1;
    }
    Sampl* uids = NULL;
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN+1];

    while (get_line_FS(sample_fp, line)) {
        if (curr == n) {
            n += 100;
            Sampl* tmp = realloc(uids, n * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_sf(): Cannot realloc.\n" );
                if (uids) {
                    free(uids);
                }
                destroy_FS(sample_fp);
                return 1;
            }
            uids = tmp;
        }
        line[strlen(line)-1] = '\0';
        Sampl* sp = find_sample(data, line);
        if (!sp) {
            fprintf( stderr, "Sample %s not found in reference panel.\n", line );
            continue;
        }
        uids[curr] = *sp;
        curr++;
    }  
    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_sf(): No matching samples found in %s.\n", sample_fn );
        free(uids);
        destroy_FS(sample_fp);
        return 1;
    }
    data->uids = uids;
    data->n_uids = curr;
    destroy_FS(sample_fp);
    return 0;
}


int read_scmd(Comp_dt* data, const char* sample_str) {
    Sampl* uids = NULL;
    int curr = 0;
    int n = 0;
    char buf[MAX_LINE_LEN+1];
    strcpy(buf, sample_str);

    char* sample = strtok(buf, ",");
    while(sample) {
        if (curr == n) {
            n += 100;
            Sampl* tmp = realloc(uids, n * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_scmd(): Cannot realloc.\n" );
                if (uids) {
                    free(uids);
                }
                return 1;
            }
            uids = tmp;
        }
        Sampl* sp = find_sample(data, sample);
        if (!sp) {
            fprintf( stderr, "Sample %s not found in reference panel.\n", sample );
            sample = strtok(NULL, ",");
            continue;
        }
        uids[curr] = *sp;
        curr++;
        sample = strtok(NULL, ",");
    }

    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_scmd(): No matching samples found.\n" );
        free(uids);
        return 1;
    }
    data->uids = uids;
    data->n_uids = curr;
    return 0;
}


int read_rf(Comp_dt* data, const char* ref_fn) {
    File_Src* ref_fp = init_FS(ref_fn);
    if (!ref_fp) {
        return 1;
    }
    Sampl* refids = NULL;
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN+1];

    while (get_line_FS(ref_fp, line)) {
        if (curr == n) {
            n += 100;
            Sampl* tmp = realloc(refids, n * sizeof(Sampl));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_rf(): Cannot realloc.\n" );
                if (refids) {
                    free(refids);
                }
                destroy_FS(ref_fp);
                return 1;
            }
            refids = tmp;
        }
        line[strlen(line)-1] = '\0';
        Sampl* sp = find_sample(data, line);
        if (!sp) {
            fprintf( stderr, "Reference sample %s not found in input panel.\n", line );
            continue;
        }
        refids[curr] = *sp;
        curr++;
    }  
    if (curr == 0) {
        fprintf( stderr, "[::] ERROR in read_rf(): No matching samples found in %s.\n", ref_fn );
        free(refids);
        destroy_FS(ref_fp);
        return 1;
    }
    data->refids = refids;
    data->n_refids = curr;
    destroy_FS(ref_fp);
    return 0;
}


int read_af(Comp_dt* data, const char* af_fn, const char* chr) {
    File_Src* af_fp = init_FS(af_fn);
    if (!af_fp) {
        return 1;
    }
    Freq* uaf = NULL;
    char buf_chr[128];
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN];

    while (get_line_FS(af_fp, line)) {
        if (curr == n) {
            n += 200;
            Freq* tmp = realloc(uaf, n * sizeof(Freq));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_af(): Cannot realloc.\n" );
                if (uaf) {
                    free(uaf);
                }
                destroy_FS(af_fp);
                return 1;
            }
            uaf = tmp;
        }
        if (sscanf(line, "%128s %lu %lf", buf_chr, &uaf[curr].pos, &uaf[curr].f) == 3) {
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
        free(uaf);
        destroy_FS(af_fp);
        return 1;
    }
    data->uaf = uaf;
    data->n_uaf = curr;
    destroy_FS(af_fp);
    return 0;
}


int read_pos(Comp_dt* data, const char* pos_fn, const char* chr) {
    File_Src* pos_fp = init_FS(pos_fn);
    if (!pos_fp) {
        return 1;
    }
    unsigned long* upos = NULL;
    char buf_chr[128];
    int curr = 0;
    int n = 0;
    char line[MAX_LINE_LEN];

    while (get_line_FS(pos_fp, line)) {
        if (curr == n) {
            n += 2000;
            unsigned long* tmp = realloc(upos, n * sizeof(unsigned long));
            if (!tmp) {
                fprintf( stderr, "[::] ERROR in read_pos(): Cannot realloc.\n" );
                if (upos) {
                    free(upos);
                }
                destroy_FS(pos_fp);
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
        destroy_FS(pos_fp);
        return 1;
    }

    data->upos = upos;
    data->n_upos = curr;
    destroy_FS(pos_fp);
    return 0;
}


int site_in_upos(Comp_dt* data, const size_t pos) {
    for (int i = 0; i < data->n_upos; i++) {
        if (pos == data->upos[i]) {
            return 1;
        }
    }
    return 0;
}


int destroy_CD(Comp_dt* data) {
    if (!data) {
        return 0;
    }
    free(data->ids);
    if (data->uids != data->ids) {
        free(data->uids);
    }
    if (data->refids != data->ids) {
        free(data->refids);
    }
    free(data->upos);
    free(data->uaf);
    free(data);
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
