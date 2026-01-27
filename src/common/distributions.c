#include "distributions.h"
#include "utils.h"
#include "structure_alias.h"

#include "vec_znx_arithmetic_private.h"


int new_uniform_random_vec(uint64_t limb_len, 
                           int64_t* res, int64_t res_size, int64_t res_sl, 
                           uint64_t nb_bits)
{
    for(uint64_t i = 0; i < res_size; i++){
        for(uint64_t j = 0; j < limb_len ; j++){
            if(rand_uniform(res + (i * res_sl) + j, nb_bits) < 0){
                return -1;}}}
    return 0;
}

int new_uniform_random_vec_dft(MODULE* module, 
                               PolyUnivDFT* res_dft, int64_t res_size,
                               uint64_t nb_bits)
{
    uint64_t N = module->nn;
    int64_t* tmp_space = malloc(N * res_size * sizeof(int64_t));
    for(int i = 0; i < res_size; i++){
        for(int p = 0; p < N ; p++){
            if(rand_uniform(tmp_space + (i*N) + p, nb_bits) < 0){
                free(tmp_space);
                return -1;}}}

    vec_znx_dft_p(module, res_dft, res_size, tmp_space, res_size, N);
    
    free(tmp_space);

    return 0;
}

int new_normal_random_vec(uint64_t limb_len, 
                      double* res, int64_t res_size, int64_t res_sl,
                      double mu, double sigma)
{
    for(int i = 0; i < res_size; i++)
        for(int j = 0; j < limb_len ; j++)
            if(rand_normal(res + (i * res_sl) + j, mu, sigma) < 0) 
                return -1;
    return 0;
}
