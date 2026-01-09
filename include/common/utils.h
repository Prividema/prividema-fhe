#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

int rand_uniform(int64_t *result);
int rand_normal(double *result, double mu, double sigma);

#endif  // UTILS_H