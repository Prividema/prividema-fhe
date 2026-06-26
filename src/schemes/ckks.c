#include "ckks.h"

#include <math.h>

void encode_slow_internal(double* out, uint64_t n, double complex* in)
{
	for (uint64_t i = 0; i < n; ++i)
	{
		double complex sum = 0;
		uint64_t pw        = 1;
		uint64_t pw_mask   = 2 * n - 1;
		for (uint64_t j = 0; j < n / 2; ++j)
		{
			double complex expo = pw * i * (M_PI * I) / n;
			sum += cexp(expo) * in[j];
			pw *= 5;
			pw &= pw_mask;
		}
		out[i] = 2.0 / n * creal(sum);
	}
}

void decode_slow_internal(double complex* out, uint64_t n, double* in)
{
	uint64_t pw      = 1;
	uint64_t pw_mask = 2 * n - 1;
	for (uint64_t j = 0; j < n / 2; ++j)
	{
		double complex expo = pw * (M_PI * I) / n;
		double complex root = cexp(-expo);
		double complex sum  = 0;
		for (int64_t i = (int64_t)n - 1; i >= 0; --i)
		{
			sum *= root;
			sum += in[i];
		}
		out[j] = sum;
		pw *= 5;
		pw &= pw_mask;
	}
}
