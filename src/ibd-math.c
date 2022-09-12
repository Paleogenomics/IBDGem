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
