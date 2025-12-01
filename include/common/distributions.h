#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H

#include <stdint.h>
#include "utils.h"

int uniform_random_vec(int64_t limb_size, int64_t* res, int64_t res_size, int64_t res_sl);

#endif // DISTRIBUTION_H