#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "load_i2.h"
#include "file-io.h"

int main(int argc, char* argv[]) {

    Impute2* data = init_I2(argv[1], argv[2], argv[3]);
    int **haps = data->haps;
    printf("Num sites: %zu\n", data->num_sites);
    printf("Num haps: %zu\n\n", data->num_haps);

    printf("Line 1 - ID: %s\n", data->ids[0]);
    printf("Line 1 - Pos: %s\n", data->pos[0]);
    printf("Line 1 - Ref: %s\n", data->ref_alleles[0]);
    printf("Line 1 - Alt: %s\n\n", data->alt_alleles[0]);

    printf("Line 100 - ID: %s\n", data->ids[99]);
    printf("Line 100 - Pos: %s\n", data->pos[99]);
    printf("Line 100 - Ref: %s\n", data->ref_alleles[99]);
    printf("Line 100 - Alt: %s\n\n", data->alt_alleles[99]);
    
    printf("Sample 1: %s\n", data->samples[0]);
    printf("Sample 100: %s\n\n", data->samples[99]);

    printf("Individual %s has alleles %d and %d at variant %s\n", data->samples[0], *(haps[0]+0), *(haps[0]+1), data->ids[0]);

    printf("Data loaded successfully.\n");
    destroy_I2(data);
}
