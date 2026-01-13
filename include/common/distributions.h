#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>
#include "utils.h"

int uniform_random_vec(int64_t limb_len, int64_t* res, int64_t n_limbs, int64_t res_sl, int64_t min, int64_t max);
int normal_random_vec(int64_t limb_len, int64_t* res, int64_t n_limbs, int64_t res_sl);

#endif // DISTRIBUTION_H