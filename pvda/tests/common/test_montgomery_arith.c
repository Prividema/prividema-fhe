#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdint.h>

#include "montgomery_arith.h"

// Tests of montgomery_arith.h (radix R = 2^32), checked against plain 128-bit modular arithmetic.

// Odd moduli below the 2^31 limit, from the smallest one to 2^31 - 1.
static const uint64_t moduli[] = {3, 5, 97, 7681, 12289, 65537, 1073741827, 2147483647};
#define NB_MODULI (sizeof(moduli) / sizeof(moduli[0]))

static const uint64_t R_BITS = 32;

/** @brief Values to try modulo q: the edges of [0, q) and a few others. */
static uint64_t sample_value(uint64_t q, int i)
{
	switch (i)
	{
		case 0:
			return 0;
		case 1:
			return 1;
		case 2:
			return q - 1;
		case 3:
			return q / 2;
		default:
			return (UINT64_C(0x9e3779b97f4a7c15) * (uint64_t)i) % q;
	}
}
#define NB_SAMPLES 40

Test(montgomery, tild_is_minus_inverse_of_q)
{
	for (size_t i = 0; i < NB_MODULI; ++i)
	{
		uint64_t q      = moduli[i];
		uint64_t q_tild = montgomery_tild_32bit(q);
		cr_assert(lt(u64, q_tild, UINT64_C(1) << R_BITS));
		// q * q_tild = -1 mod 2^32
		cr_assert(eq(u64, (q * q_tild) & 0xFFFFFFFFu, 0xFFFFFFFFu), "q = %lu", q);
	}
}

Test(montgomery, r2modq)
{
	for (size_t i = 0; i < NB_MODULI; ++i)
	{
		uint64_t q        = moduli[i];
		unsigned __int128 r2 = (unsigned __int128)1 << 64;
		cr_assert(eq(u64, montgomery_r2modq_32bit(q), (uint64_t)(r2 % q)), "q = %lu", q);
	}
}

Test(montgomery, red_is_division_by_R)
{
	for (size_t i = 0; i < NB_MODULI; ++i)
	{
		uint64_t q      = moduli[i];
		uint64_t q_tild = montgomery_tild_32bit(q);
		// Inputs in [0, q * R)
		const uint64_t inputs[] = {0, 1, q - 1, q, (q << R_BITS) - 1, (q << R_BITS) / 3, UINT64_C(12345678901) % (q << R_BITS)};
		for (size_t j = 0; j < sizeof(inputs) / sizeof(inputs[0]); ++j)
		{
			uint64_t r = montgomery_red_32bit(inputs[j], q, q_tild);
			cr_assert(lt(u64, r, q), "q = %lu, m = %lu", q, inputs[j]);
			// r * R = m mod q
			cr_assert(eq(u64, (uint64_t)(((unsigned __int128)r << R_BITS) % q), inputs[j] % q), "q = %lu, m = %lu",
			          q, inputs[j]);
		}
	}
}

Test(montgomery, encode_decode_round_trip)
{
	for (size_t i = 0; i < NB_MODULI; ++i)
	{
		uint64_t q      = moduli[i];
		uint64_t q_tild = montgomery_tild_32bit(q);
		uint64_t r2     = montgomery_r2modq_32bit(q);
		for (int s = 0; s < NB_SAMPLES; ++s)
		{
			uint64_t x   = sample_value(q, s);
			uint64_t x_m = montgomery_encode_32bit(x, q, q_tild, r2);
			cr_assert(lt(u64, x_m, q));
			// x_m = x * R mod q
			cr_assert(eq(u64, x_m, (uint64_t)(((unsigned __int128)x << R_BITS) % q)), "q = %lu, x = %lu", q, x);
			cr_assert(eq(u64, montgomery_decode_32bit(x_m, q, q_tild), x), "q = %lu, x = %lu", q, x);
		}
	}
}

Test(montgomery, mult_matches_modular_product)
{
	for (size_t i = 0; i < NB_MODULI; ++i)
	{
		uint64_t q      = moduli[i];
		uint64_t q_tild = montgomery_tild_32bit(q);
		uint64_t r2     = montgomery_r2modq_32bit(q);
		for (int a = 0; a < NB_SAMPLES; ++a)
			for (int b = 0; b < 8; ++b)
			{
				uint64_t x = sample_value(q, a), y = sample_value(q, b);
				uint64_t p_m = montgomery_mult_32bit(montgomery_encode_32bit(x, q, q_tild, r2),
				                                     montgomery_encode_32bit(y, q, q_tild, r2), q, q_tild);
				cr_assert(eq(u64, montgomery_decode_32bit(p_m, q, q_tild), (x * y) % q), "q = %lu, %lu * %lu", q, x,
				          y);
			}
	}
}

Test(montgomery, pow_exp)
{
	for (size_t i = 0; i < NB_MODULI; ++i)
	{
		uint64_t q      = moduli[i];
		uint64_t q_tild = montgomery_tild_32bit(q);
		uint64_t r2     = montgomery_r2modq_32bit(q);
		uint64_t one_m  = montgomery_encode_32bit(1, q, q_tild, r2);
		for (int s = 0; s < 10; ++s)
		{
			uint64_t x   = sample_value(q, s);
			uint64_t x_m = montgomery_encode_32bit(x, q, q_tild, r2);

			// x^0 = 1 (also for x = 0), x^1 = x
			cr_assert(eq(u64, montgomery_decode_32bit(montgomery_pow_exp_32bit(x_m, 0, q, q_tild, one_m), q, q_tild), 1 % q));
			cr_assert(eq(u64, montgomery_decode_32bit(montgomery_pow_exp_32bit(x_m, 1, q, q_tild, one_m), q, q_tild), x));

			// Against naive exponentiation
			const uint64_t exps[] = {2, 3, 17, 1000, q - 1};
			for (size_t e = 0; e < sizeof(exps) / sizeof(exps[0]); ++e)
			{
				uint64_t expected = 1 % q;
				for (uint64_t b = x, n = exps[e]; n; n >>= 1, b = (b * b) % q)
					if (n & 1) expected = (expected * b) % q;
				uint64_t got =
				    montgomery_decode_32bit(montgomery_pow_exp_32bit(x_m, exps[e], q, q_tild, one_m), q, q_tild);
				cr_assert(eq(u64, got, expected), "q = %lu, %lu ^ %lu", q, x, exps[e]);
			}
		}
	}
}
