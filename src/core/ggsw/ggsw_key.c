#include "ggsw_key.h"
#include "stdio.h"

GGSWSecretKeyDFT* new_ggsw_secret_key_dft(PolyUnivDFT** values, uint64_t N, uint64_t k
){
    GGSWSecretKeyDFT* sk_dft = malloc(sizeof(GGSWSecretKeyDFT));
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

    return sk_dft;
}

void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft
){
    
    for(int64_t j = 0 ; j < sk_dft->k ; j++)
    {
        free(sk_dft->values[j]);
    }
    free(sk_dft->values);
    free(sk_dft);
}