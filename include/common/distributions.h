#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>
#include "utils.h"

int uniform_random_vec(int64_t* res, int64_t res_size, gsl_rng *r);

#endif // DISTRIBUTION_H