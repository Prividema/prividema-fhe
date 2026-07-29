#include "math_utils.h"

#include <math.h>
#include <stdint.h>

uint64_t mod_inv(int64_t x, int64_t mod)
{
	// Extended Euclidean algorithm
	// See https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm#Computing_multiplicative_inverses_in_modular_structures
	// if you need a refresher/source
	int64_t r, t, r_new, t_new;
	t     = 0;
	r     = mod;
	t_new = 1;
	r_new = x;

	while (r_new)
	{
		uint64_t q = r / r_new;

		uint64_t t_tmp = t - q * t_new;
		t              = t_new;
		t_new          = t_tmp;

		uint64_t r_tmp = r - q * r_new;
		r              = r_new;
		r_new          = r_tmp;
	}

	if (t < 0) t += mod;
	return t;
}

double rnx_torus_distance(double a, double b)
{
	a -= floor(a);
	b -= floor(b);
	if (a > b)
	{
		return fmin(a - b, (b + 1) - a);
	}
	return fmin(b - a, (a + 1) - b);
}

uint64_t tnx_torus_distance(uint64_t a, uint64_t b) { return a - b > b - a ? b - a : a - b; }

uint64_t u64_round_up_div(uint64_t num, uint64_t div) { return (num + (div - 1)) / div; }
int32_t i32_round_up_div(int32_t num, int32_t div) { return (num + (div - 1)) / div; }
uint64_t next_pow2_log(uint64_t num) { return num <= 1 ? 0 : 64 - __builtin_clzll(num - 1); }
