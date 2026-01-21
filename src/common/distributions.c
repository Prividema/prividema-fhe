#include "distributions.h"
#include "utils.h"

int uniform_random_vec(int64_t limb_len, 
                       int64_t* res, int64_t res_size, int64_t res_sl, 
                       int nb_bits)
{
    for(int i = 0; i < res_size; i++)
        for(int j = 0; j < limb_len ; j++)
            if(rand_uniform(res + (i * res_sl) + j, nb_bits) < 0)
                return -1;
    return 0;
}


int normal_random_vec(int64_t limb_len, 
                      double* res, int64_t res_size, int64_t res_sl,
                      double mu, double sigma)
{
    for(int i = 0; i < res_size; i++)
        for(int j = 0; j < limb_len ; j++)
            if(rand_normal(res + (i * res_sl) + j, mu, sigma) < 0) 
                return -1;
    return 0;
}
