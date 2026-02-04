#include "glwe_key.h"
#include "rng.h"
#include "spqlios_alias.h"

//! GLWE PART (begin)

PolyUniv** new_glwe_secret_key_values(uint64_t N, uint64_t k
){
    PolyUniv** values = malloc(k*sizeof(PolyUniv*));

    for(int64_t j = 0 ; j < k ; j++)
        values[j] = calloc(N, sizeof(int64_t));
    
    return values;
}   

PolyUniv** new_uniform_glwe_secret_key_values(uint64_t N, uint64_t k, uint64_t  nb_bits
){
    PolyUniv** values = malloc(k*sizeof(PolyUniv*));

    // Uniform random generation of k Zn[X] polynomials.
    MODULE* module = new_module_info_p(N);
    for(int64_t j = 0 ; j < k ; j++){
        values[j] = calloc(N, sizeof(int64_t));
        inplace_uniform_random_vec(N, values[j], 1, N, nb_bits);
    }
    
    delete_module_info_p(module);

    return values;
}   

void delete_glwe_secret_key_values(PolyUniv**values, uint64_t k){
    for(int64_t j = 0 ; j < k ; j++){
        free(values[j]);}
    
    free(values);
}

GLWESecretKey* new_glwe_secret_key(PolyUniv** values, uint64_t N, uint64_t k
){
    GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
    if(sk == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk->N = N;
    sk->k = k;

    if(values == NULL)
    {
        sk->values = new_glwe_secret_key_values(N, k);
    }
    else
    {
        sk->values = values;
    }

    return sk;
}

GLWESecretKey* new_uniform_glwe_secret_key(uint64_t N, uint64_t k, uint64_t nb_bits){
    GLWESecretKey* sk = malloc(sizeof(GLWESecretKey));
    if(sk == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk->N = N;
    sk->k = k;
    sk->values = new_uniform_glwe_secret_key_values(N, k, nb_bits);

    return sk;
}

void delete_glwe_secret_key(GLWESecretKey* sk
){
    delete_glwe_secret_key_values(sk->values, sk->k);
    free(sk);
}

GLWESecretKey* transform_glwe_secret_key_dft_to_not_dft(GLWESecretKeyDFT* sk_dft
){
    GLWESecretKey* sk = new_glwe_secret_key(NULL, sk_dft->N, sk_dft->k);
    MODULE* module = new_module_info_p(sk_dft->N);

    for(int64_t j = 0 ; j < sk->k ; j++){
        vec_znx_idft_p(module, sk->values[j], 1, sk_dft->values[j], 1);
    }

    delete_module_info_p(module);

    return sk;
}

PolyUniv** transform_glwe_secret_key_values_dft_to_not_dft(PolyUnivDFT** values_dft, uint64_t N, uint64_t k
){
    PolyUniv** values = malloc(k*sizeof(PolyUniv*));
    MODULE* module = new_module_info_p(N);

    for(int64_t j = 0 ; j < k ; j++)
    {
        values[j] = calloc(N, sizeof(int64_t));
        vec_znx_idft_p(module, values[j], 1, values_dft[j], 1);
    }
    
    delete_module_info_p(module);

    return values;
}

//! GLWE PART in DFT space (begin)

PolyUnivDFT** new_glwe_secret_key_values_dft(uint64_t N, uint64_t k
){
    PolyUnivDFT** values = malloc(k*sizeof(PolyUnivDFT*));

    for(int64_t j = 0 ; j < k ; j++){
        values[j] = calloc(N, sizeof(double));
    }

    return values;
}   

PolyUnivDFT** new_uniform_glwe_secret_key_values_dft(uint64_t N, uint64_t k, uint64_t nb_bits
){
    PolyUnivDFT** values_dft = malloc(k*sizeof(PolyUnivDFT*));

    // Uniform random generation of k Zn[X] polynomials.
    MODULE* module = new_module_info_p(N);
    for(int64_t j = 0 ; j < k ; j++){
        values_dft[j] = malloc(N * sizeof(double));
        if(inplace_uniform_random_vec_znx_dft(module, values_dft[j], 1, nb_bits) < 0){
            delete_module_info_p(module);
            return NULL;
        }
        
    }
    
    delete_module_info_p(module);
    
    return values_dft;
}   

void delete_glwe_secret_key_values_dft(PolyUnivDFT**values, uint64_t k){
    for(int64_t j = 0 ; j < k ; j++)
        free(values[j]);
    
    free(values);
}

GLWESecretKeyDFT* new_glwe_secret_key_dft(PolyUnivDFT** values, uint64_t N, uint64_t k
){
    GLWESecretKeyDFT* sk_dft = malloc(sizeof(GLWESecretKeyDFT));
    if(sk_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk_dft->N = N;
    sk_dft->k = k;

    if(values == NULL)
    {
        sk_dft->values = new_glwe_secret_key_values_dft(N, k);
    }
    else
    {
        sk_dft->values = values;
    }

    return sk_dft;
}

GLWESecretKeyDFT* new_uniform_glwe_secret_key_dft(uint64_t N, uint64_t k, uint64_t nb_bits
){
    GLWESecretKeyDFT* sk_dft = malloc(sizeof(GLWESecretKeyDFT));
    if(sk_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk_dft->N = N;
    sk_dft->k = k;
    sk_dft->values = new_uniform_glwe_secret_key_values_dft(N, k, nb_bits);

    return sk_dft;
}

void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft
){
    delete_glwe_secret_key_values_dft(sk_dft->values, sk_dft->k);
    free(sk_dft);
}

GLWESecretKeyDFT* transform_glwe_secret_key_not_dft_to_dft(GLWESecretKey* sk
){
    GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(NULL, sk->N, sk->k);
    MODULE* module = new_module_info_p(sk_dft->N);

    for(int64_t j = 0 ; j < sk_dft->k ; j++){
        vec_znx_dft_p(module, sk_dft->values[j], 1, sk->values[j], 1, sk->N);
    }

    delete_module_info_p(module);

    return sk_dft;
}

PolyUnivDFT** transform_glwe_secret_key_values_not_dft_to_dft(PolyUniv** values, uint64_t N, uint64_t k
){
    PolyUnivDFT** values_dft = malloc(k*sizeof(PolyUniv*));
    MODULE* module = new_module_info_p(N);

    for(int64_t j = 0 ; j < k ; j++)
    {
        values_dft[j] = calloc(N, sizeof(double));
        vec_znx_dft_p(module, values_dft[j], 1, values[j], 1, N);
    }
    
    delete_module_info_p(module);

    return values_dft;
}