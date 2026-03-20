#include "utils.h"

#include <math.h>
#include <stdint.h>

double binomial_tail(uint64_t N, double p, int k)
{
	double mu    = N * p;
	double sigma = sqrt(N * p * (1.0 - p));

	double r = ceil(mu * k);

	double z = (r - 0.5 - mu) / sigma;

	return 0.5 * erfc(z / sqrt(2.0));
}

double torus_distance(double a, double b)
{
	a -= floor(a);
	b -= floor(b);
	if (a > b)
	{
		return fmin(a - b, (b + 1) - a);
	}
	return fmin(b - a, (a + 1) - b);
}
