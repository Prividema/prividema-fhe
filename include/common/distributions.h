#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>
#include "utils.h"

int uniform_random_vec(int64_t limb_len, int64_t* res, int64_t res_size, int64_t res_sl, int nb_bits);
int normal_random_vec(int64_t limb_len, double* res, int64_t res_size, int64_t res_sl, double mu, double sigma);

#endif // DISTRIBUTION_H