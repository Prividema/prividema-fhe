#include "glwe_key.h"
#include <stdio.h>

GLWEPreparedSK* new_glwe_secret_key_dft(PolyUnivDFT** values, int64_t N, int64_t k){
    GLWEPreparedSK* sk_dft = malloc(sizeof(GLWEPreparedSK));
    if(sk_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk_dft->N = N;
    sk_dft->k = k;

    if(values == NULL)
    {
        sk_dft->values = malloc(k*sizeof(PolyUnivDFT*));
        for(int64_t j = 0 ; j < k ; j++)
        {
            sk_dft->values[j] = malloc(N * sizeof(int64_t));
        }
    }
    else
    {
        sk_dft->values = values;
    }
}