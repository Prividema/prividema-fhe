#include "heint.h"

#include <assert.h>
#include <stdint.h>
#include <sys/types.h>

const uint64_t bitmask_32bit = 0x00000000FFFFFFFF;

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

uint64_t montgomery_encode_32bit(uint64_t x, uint64_t q) { return (x << 32) % q; }
uint64_t montgomery_decode_32bit(uint64_t x, uint64_t q, uint64_t q_tild) { return montgomery_red_32bit(x, q, q_tild); }
uint64_t montgomery_mult_32bit(uint64_t x_m, uint64_t y_m, uint64_t q, uint64_t q_tild)
{
	return montgomery_red_32bit(x_m * y_m, q, q_tild);
}

uint64_t montgomery_red_32bit(uint64_t m, uint64_t q, uint64_t q_tild)
{
	assert(q <= (1ll << 31));
	uint64_t v = ((m & bitmask_32bit) * q_tild) & bitmask_32bit;
	uint64_t w = m + v * q;
	uint64_t r = w >> 32;
	/*
   * The following is a branchless version of:
	  if (r >= q) return r - q;
	  return r;
	*/
	uint64_t msk = -((uint64_t)(r >= q));  // All 0 or all 1
	return r - (msk & q);
}

uint64_t montgomery_pow_exp_32bit(uint64_t base_m, uint64_t exp, uint64_t q, uint64_t q_tild, uint64_t one_m)
{
	uint64_t ans_m = one_m;
	uint64_t mult  = base_m;

	while (exp > 0)
	{
		if (exp & 1)
		{
			ans_m = montgomery_mult_32bit(ans_m, mult, q, q_tild);
		}
		exp  = exp / 2;
		mult = montgomery_mult_32bit(mult, mult, q, q_tild);
	}

	return ans_m;
}
