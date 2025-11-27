#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_rng.h>

void gsl_init(gsl_rng *r);
int gsl_refresh_seed(gsl_rng *r);

#endif  // UTILS_H