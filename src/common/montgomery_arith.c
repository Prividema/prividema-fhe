#include "montgomery_arith.h"

#include <assert.h>

#include "math_utils.h"

const uint64_t bitmask_32bit = 0x00000000FFFFFFFF;

uint64_t montgomery_encode_32bit(uint64_t x, uint64_t q, uint64_t q_tild, uint64_t r2modn)
{
	// Constant time encoding by REDC(x * ((r^2) mod q))
	return montgomery_mult_32bit(x, r2modn, q, q_tild);
}
uint64_t montgomery_tild_32bit(uint64_t q) { return (1ull << 32) - mod_inv(q, 1ull << 32); }

uint64_t montgomery_r2modq_32bit(uint64_t q)
{
	// We want to compute r^2 mod q, but that would overflow a single value
	// since r = 2^32 ==> r^2 = 2^64
	// Therefore we want 2^64 mod q, which is equal to (2^64 - q) mod q,
	// and which in C unsigned arithmetic is represented by (-q) mod q
	return (-q) % q;
}
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

	// square and multiply
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
