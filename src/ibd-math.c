#include "ibd-math.h"


/* Private function for performing the nCk calculation */
static unsigned long choose(unsigned int n, unsigned int k) {
    if (k == 0) {
        return 1;
    }
    return (n * choose(n-1, k-1)) / k;
}


unsigned long** init_nCk(unsigned int n) {
    unsigned long** nCk = malloc((n+1) * sizeof(unsigned long*));
    for (int i = 0; i <= n; i++) {
        unsigned long* iCk = malloc((n+1) * sizeof(unsigned long));
        for (int j = 0; j <= n; j++) {
            iCk[j] = choose(i, j);
        }
        nCk[i] = iCk;
    }  
    return nCk;
}


unsigned long retrieve_nCk(unsigned long** nCk, unsigned int n, unsigned int k) {
    if (k > n) {
        fprintf(stderr, "Invalid k = %u; enter k so that n>=k>=0\n", k);
    }
    return *(nCk[n]+k);
}


int destroy_nCk(unsigned long** nCk, unsigned int n) {
    if (!nCk) {
        return 0;
    }
    for (int i = 0; i <= n; i++) {
        free(nCk[i]);
    }
    free(nCk);
    return 0;
}


double find_pDgG(unsigned long** nCk, double epsilon, unsigned short A0, unsigned short A1,
                 unsigned int n_ref, unsigned int n_alt) {
    double p = 1.0; 
    
    // probability is 1 if there are no observed bases
    if (n_ref == 0 && n_alt == 0) {
        return p;
    }

    unsigned long nCk_val = retrieve_nCk(nCk, n_ref+n_alt, n_ref);
    
    if (A0 == 0 && A1 == 0) {
        p = nCk_val * pow(1-epsilon, n_ref) * pow(epsilon, n_alt);
    }
    else if (A0 == 1 && A1 == 1) {
        p = nCk_val * pow(1-epsilon, n_alt) * pow(epsilon, n_ref);
    }

    // If heterozygous site, prob of drawing either allele is 0.5
    // -> binomial function with p=0.5
    // (n choose k) x p**k x (1-p)**(n-k)
    // here n_ref is arbitrarily chosen to be 'successes'
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        p = nCk_val * pow(0.5, n_ref) * pow(0.5, n_alt);
    }
    else {
        fprintf( stderr, "Invalid genotype: A0 = %d, A1 = %d\n", A0, A1 );
        exit(1);
    }

    // avoids multiplying with 0 during likelihood aggregation
    if (p == 0.0) {
        p = DBL_MIN; 
    }
    return p;
}


double find_pDgf(double f, double pD_g_00, double pD_g_01, double pD_g_11) {
    double p = 1.0;
    
    // only case where these probs are 1.0 is when nREF+nALT = 0
    if (pD_g_00 == 1 || pD_g_01 == 1 || pD_g_11 == 1) {
        return p;
    }

    // Hardy-Weinberg  P(Data | Genotype)
    p = pow(1-f, 2.0)  *  pD_g_00 +
        2 * (1-f) * f  *  pD_g_01 +
        pow(f, 2.0)    *  pD_g_11;

    if (p == 0.0) {
        p = DBL_MIN;
    }
    return p;
}


double find_pDgIBD1(unsigned short A0, unsigned short A1, double f,
                    double pD_g_00, double pD_g_01, double pD_g_11) {
    double p = 1.0;

    // For each VCF genotype, use the probability of the *underlying genotypes* 
    // of the BAM data to calculate the probability of the BAM data given that
    // VCF genotype under IBD1 model

    // The probabilities of all 3 underlying BAM genotypes ((0,0), (0,1), (1,1))
    // can be found for any given VCF genotype 
    
    // Homozygous REF: P(G_BAM=(0,1) | G_VCF=(0,0)) * P(D | G_BAM=(0,1)) +
    //                 P(G_BAM=(0,0) | G_VCF=(0,0)) * P(D | G_BAM=(0,0)) + 0
    // -> 0 prob of seeing BAM genotype (1,1) given VCF genotype (0,0)
    if (A0 == 0 && A1 == 0) {
        p = (f * pD_g_01) + ((1-f) * pD_g_00);
    }

    // Heterozygous: P(G_BAM=(0,1) | G_VCF=(0,1)) * P(D | G_BAM=(0,1)) +
    //               P(G_BAM=(0,0) | G_VCF=(0,1)) * P(D | G_BAM=(0,0)) +
    //               P(G_BAM=(1,1) | G_VCF=(0,1)) * P(D | G_BAM=(1,1))
    else if ((A0 == 0 && A1 == 1) || (A0 == 1 && A1 == 0)) {
        p = (0.5 *         pD_g_01) +
            (0.5 * (1-f) * pD_g_00) +
            (0.5 *   f   * pD_g_11);
    }

    // Homozygous ALT: P(G_BAM=(0,1) | G_VCF=(1,1)) * P(D | G_BAM=(0,1)) +
    //                 P(G_BAM=(1,1) | G_VCF=(1,1)) * P(D | G_BAM=(1,1)) + 0
    // -> 0 prob of seeing BAM genotype (0,0) given VCF genotype (1,1)
    else if (A0 == 1 && A1 == 1) {
        p = ((1-f) * pD_g_01) + (f * pD_g_11);
    }

    if (p == 0.0) {
        p = DBL_MIN;
    }
    return p;
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