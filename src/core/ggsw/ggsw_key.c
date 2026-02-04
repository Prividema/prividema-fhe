#include "ggsw_key.h"
#include "rng.h"
#include "spqlios_alias.h"

//! GGSW PART (begin)

PolyUniv** new_ggsw_secret_key_values(uint64_t N, uint64_t k
){
    PolyUniv** values = malloc(k*sizeof(PolyUniv*));

    for(int64_t j = 0 ; j < k ; j++)
        values[j] = calloc(N, sizeof(int64_t));
    
    return values;
}   

PolyUniv** new_uniform_ggsw_secret_key_values(uint64_t N, uint64_t k, uint64_t  nb_bits
){
    PolyUniv** values = malloc(k*sizeof(PolyUniv*));

    // Uniform random generation of k Zn[X] polynomials.
    MODULE* module = new_module_info_p(N);
    for(int64_t j = 0 ; j < k ; j++){
        values[j] = malloc(N * sizeof(double));
        inplace_uniform_random_vec(N, values[j], 1, N, nb_bits);
    }
    
    delete_module_info_p(module);
    
    return values;
}   

void delete_ggsw_secret_key_values(PolyUniv**values, uint64_t k){
    for(int64_t j = 0 ; j < k ; j++)
        free(values[j]);
    
    free(values);
}

GGSWSecretKey* new_ggsw_secret_key(PolyUniv** values, uint64_t N, uint64_t k
){
    GGSWSecretKey* sk = malloc(sizeof(GGSWSecretKey));
    if(sk == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk->N = N;
    sk->k = k;

    if(values == NULL)
    {
        sk->values = new_ggsw_secret_key_values(N, k);
    }
    else
    {
        sk->values = values;
    }

    return sk;
}

GGSWSecretKey* new_uniform_ggsw_secret_key(uint64_t N, uint64_t k, uint64_t nb_bits){
    GGSWSecretKey* sk = malloc(sizeof(GGSWSecretKey));
    if(sk == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk->N = N;
    sk->k = k;
    sk->values = new_uniform_ggsw_secret_key_values(N, k, nb_bits);

    return sk;
}

void delete_ggsw_secret_key(GGSWSecretKey* sk
){
    delete_ggsw_secret_key_values(sk->values, sk->k);
    free(sk);
}

GGSWSecretKey* transform_ggsw_secret_key_dft_to_not_dft(GGSWSecretKeyDFT* sk_dft
){
    GGSWSecretKey* sk = new_ggsw_secret_key(NULL, sk_dft->N, sk_dft->k);
    MODULE* module = new_module_info_p(sk_dft->N);

    for(int64_t j = 0 ; j < sk->k ; j++){
        vec_znx_idft_p(module, sk->values[j], 1, sk_dft->values[j], 1);
    }

    delete_module_info_p(module);

    return sk;
}

PolyUniv** transform_ggsw_secret_key_values_dft_to_not_dft(PolyUnivDFT** values_dft, uint64_t N, uint64_t k
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

//! GGSW PART in DFT space (begin)

PolyUnivDFT** new_ggsw_secret_key_values_dft(uint64_t N, uint64_t k
){
    PolyUnivDFT** values = malloc(k*sizeof(PolyUnivDFT*));

    for(int64_t j = 0 ; j < k ; j++){
        values[j] = calloc(N, sizeof(double));
    }

    return values;
}   

PolyUnivDFT** new_uniform_ggsw_secret_key_values_dft(uint64_t N, uint64_t k, uint64_t nb_bits
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

void delete_ggsw_secret_key_values_dft(PolyUnivDFT**values, uint64_t k){
    for(int64_t j = 0 ; j < k ; j++)
        free(values[j]);
    
    free(values);
}

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
        sk_dft->values = new_ggsw_secret_key_values_dft(N, k);
    }
    else
    {
        sk_dft->values = values;
    }

    return sk_dft;
}

GGSWSecretKeyDFT* new_uniform_ggsw_secret_key_dft(uint64_t N, uint64_t k, uint64_t nb_bits
){
    GGSWSecretKeyDFT* sk_dft = malloc(sizeof(GGSWSecretKeyDFT));
    if(sk_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    sk_dft->N = N;
    sk_dft->k = k;
    sk_dft->values = new_uniform_ggsw_secret_key_values_dft(N, k, nb_bits);

    return sk_dft;
}

void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft
){
    delete_ggsw_secret_key_values_dft(sk_dft->values, sk_dft->k);
    free(sk_dft);
}

GGSWSecretKeyDFT* transform_ggsw_secret_key_not_dft_to_dft(GGSWSecretKey* sk
){
    GGSWSecretKeyDFT* sk_dft = new_ggsw_secret_key_dft(NULL, sk->N, sk->k);
    MODULE* module = new_module_info_p(sk_dft->N);

    for(int64_t j = 0 ; j < sk_dft->k ; j++){
        vec_znx_dft_p(module, sk_dft->values[j], 1, sk->values[j], 1, sk->N);
    }

    delete_module_info_p(module);

    return sk_dft;
}

PolyUnivDFT** transform_ggsw_secret_key_values_not_dft_to_dft(PolyUniv** values, uint64_t N, uint64_t k
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