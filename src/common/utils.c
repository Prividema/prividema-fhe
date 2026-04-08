#include "utils.h"

#include <math.h>
#include <stdint.h>

double binomial_tail(uint64_t nn, double p, int k)
{
	double mu    = nn * p;
	double sigma = sqrt(nn * p * (1.0 - p));

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

uint64_t u64_round_up_div(uint64_t num, uint64_t div) { return (num + (div - 1)) / div; }
