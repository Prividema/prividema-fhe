#include "utils.h"
#include "math.h"
#include "stdio.h"

#include <math.h>
#include <stdint.h>

double binomial_tail(uint64_t N, double p)
{
    double mu = N * p;
    double sigma = sqrt(N * p * (1.0 - p));

    double k = ceil(mu);

    double z = (k - 0.5 - mu) / sigma;

    return 0.5 * erfc(z / sqrt(2.0));
}