#include "glwe_key.h"
#include "distributions.h"
#include "spqlios_alias.h"
#include <stdio.h>

//! GLWE PART (begin)
GLWESecretKey* new_glwe_secret_key(PolyUniv** values, int64_t N, int64_t k
){
    GLWESecretKey* sk = malloc(sizeof(GLWESecretKey*));
    if(sk == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk->N = N;
    sk->k = k;

    if(values == NULL)
    {
        sk->values = malloc(k*sizeof(PolyUniv*));
        for(int64_t j = 0 ; j < k ; j++)
        {
            sk->values[j] = calloc(N * sizeof(int64_t),1);
        }
    }
    else
    {
        sk->values = values;
    }

    return sk;
}

void delete_glwe_secret_key(GLWESecretKey* sk
){
    for(int64_t j = 0 ; j < sk->k ; j++)
        free(sk->values[j]);
    
    free(sk->values);
    free(sk);
}

GLWESecretKey* glwe_secret_key_dft_to_not_dft(GLWEPreparedSK* sk_dft
){
    GLWESecretKey* sk = new_glwe_secret_key(NULL, sk_dft->N, sk_dft->k);
    MODULE* module = new_module_info_p(sk->N);

    for(int64_t j = 0 ; j < sk->k ; j++){
        vec_znx_idft_p(module, sk->values[j], 1, sk_dft->values[j], 1);
    }

    delete_module_info_p(module);
    
    return sk;
}

PolyUniv** secret_key_values_dft_to_not_dft(GLWEPreparedSK* sk_dft
){
    PolyUniv** values = malloc(sk_dft->k*sizeof(PolyUniv*));
    MODULE* module = new_module_info_p(sk_dft->N);

    for(int64_t j = 0 ; j < sk_dft->k ; j++)
    {
        values[j] = calloc(sk_dft->N * sizeof(int64_t),1);
        vec_znx_idft_p(module, values[j], 1, sk_dft->values[j], 1);
    }
    
    delete_module_info_p(module);

    return values;
}

//! GLWE PART in DFT space (begin)

GLWEPreparedSK* new_glwe_secret_key_dft(PolyUnivDFT** values, int64_t N, int64_t k
){
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
            sk_dft->values[j] = calloc(N * sizeof(double),1);
        }
    }
    else
    {
        sk_dft->values = values;
    }

    return sk_dft;
}

void delete_glwe_secret_key_dft(GLWEPreparedSK* sk_dft
){
    for(int64_t j = 0 ; j < sk_dft->k ; j++)
        free(sk_dft->values[j]);
    
    free(sk_dft->values);
    free(sk_dft);
}

GLWEPreparedSK* new_uniform_glwe_secret_key_gen(int64_t N, int64_t k, int nb_bits
){
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(NULL, N, k); 

    // Uniform random generation of k Zn[X] polynomials.
    MODULE* module = new_module_info_p(N);
    for(int64_t j = 0 ; j < k ; j++){
        new_uniform_random_vec_dft(module, sk_dft->values[j], 1, nb_bits);
    }
    
    delete_module_info_p(module);

    return sk_dft;
}

GLWEPreparedSK* glwe_secret_key_not_dft_to_dft(GLWESecretKey* sk
){
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(NULL, sk->N, sk->k);
    MODULE* module = new_module_info_p(sk_dft->N);

    for(int64_t j = 0 ; j < sk_dft->k ; j++){
        vec_znx_dft_p(module, sk_dft->values[j], 1, sk->values[j], 1, sk->N);
    }

    delete_module_info_p(module);

    return sk_dft;
}