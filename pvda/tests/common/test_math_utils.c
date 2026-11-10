#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <math.h>
#include <stdint.h>

#include "math_utils.h"

// Tests of math_utils.h: modular inverse, torus distances, rounded-up division and ceil(log2).

// =============================================
// |                                           |
// |                 mod_inv                   |
// |                                           |
// =============================================

/** @brief (x * y) mod m with x possibly negative, computed on 128 bits. */
static uint64_t mul_mod(int64_t x, uint64_t y, uint64_t m)
{
	__int128 r = ((__int128)x * (__int128)y) % (__int128)m;
	if (r < 0) r += m;
	return (uint64_t)r;
}

Test(mod_inv, small_primes_every_residue)
{
	static const int64_t primes[] = {2, 3, 5, 7, 13, 257};
	for (size_t i = 0; i < sizeof(primes) / sizeof(primes[0]); ++i)
		for (int64_t x = 1; x < primes[i]; ++x)
		{
			uint64_t inv = mod_inv(x, primes[i]);
			cr_assert(lt(u64, inv, (uint64_t)primes[i]));
			cr_assert(eq(u64, mul_mod(x, inv, primes[i]), 1), "x = %ld mod %ld", x, primes[i]);
		}
}

Test(mod_inv, composite_modulus_and_large_values)
{
	// Coprime pairs, including the 32-bit Montgomery radix and a 62-bit modulus
	cr_assert(eq(u64, mul_mod(7, mod_inv(7, 1 << 20), 1 << 20), 1));
	cr_assert(eq(u64, mul_mod(2147483647, mod_inv(2147483647, INT64_C(1) << 32), INT64_C(1) << 32), 1));
	const int64_t big = (INT64_C(1) << 61) - 1;  // Mersenne prime
	cr_assert(eq(u64, mul_mod(123456789, mod_inv(123456789, big), big), 1));
	// x = 1 and x = mod - 1 are their own inverses
	cr_assert(eq(u64, mod_inv(1, 97), 1));
	cr_assert(eq(u64, mod_inv(96, 97), 96));
}

Test(mod_inv, input_outside_of_the_ring)
{
	// x is reduced modulo mod first: negative and larger-than-mod inputs give the inverse of their residue
	cr_assert(eq(u64, mod_inv(-1, 7), 6));
	cr_assert(eq(u64, mod_inv(-3, 7), mod_inv(4, 7)));
	cr_assert(eq(u64, mod_inv(8, 7), 1));
	cr_assert(eq(u64, mod_inv(7 * 1000 + 3, 7), mod_inv(3, 7)));
	cr_assert(eq(u64, mod_inv(-(7 * 1000 + 3), 7), mod_inv(4, 7)));
}

// =============================================
// |                                           |
// |             Torus distances               |
// |                                           |
// =============================================

Test(rnx_torus_distance, basic_and_wrap_around)
{
	cr_assert(eq(dbl, rnx_torus_distance(0.25, 0.25), 0.0));
	cr_assert(eq(dbl, rnx_torus_distance(0.1, 0.3), rnx_torus_distance(0.3, 0.1)));
	cr_assert(fabs(rnx_torus_distance(0.1, 0.3) - 0.2) < 1e-15);
	// Shortest way goes through 0 = 1
	cr_assert(fabs(rnx_torus_distance(0.05, 0.95) - 0.1) < 1e-15);
	// Opposite points are at the maximal distance 1/2
	cr_assert(eq(dbl, rnx_torus_distance(0.0, 0.5), 0.5));
	cr_assert(eq(dbl, rnx_torus_distance(0.75, 0.25), 0.5));
}

Test(rnx_torus_distance, integer_parts_are_ignored)
{
	cr_assert(eq(dbl, rnx_torus_distance(3.0, -2.0), 0.0));
	cr_assert(fabs(rnx_torus_distance(-0.25, 0.75)) < 1e-15);
	cr_assert(fabs(rnx_torus_distance(10.1, -3.1) - 0.2) < 1e-12);
}

Test(tnx_torus_distance, basic_and_wrap_around)
{
	cr_assert(eq(u64, tnx_torus_distance(0, 0), 0));
	cr_assert(eq(u64, tnx_torus_distance(5, 2), 3));
	cr_assert(eq(u64, tnx_torus_distance(2, 5), 3));
	// Through 0: 2^64 - 1 and 1 are 2 apart
	cr_assert(eq(u64, tnx_torus_distance(UINT64_MAX, 1), 2));
	cr_assert(eq(u64, tnx_torus_distance(1, UINT64_MAX), 2));
	// Opposite points: distance 2^63 both ways
	cr_assert(eq(u64, tnx_torus_distance(0, UINT64_C(1) << 63), UINT64_C(1) << 63));
	cr_assert(eq(u64, tnx_torus_distance(UINT64_C(1) << 62, UINT64_C(3) << 62), UINT64_C(1) << 63));
}

// =============================================
// |                                           |
// |            Rounded-up division            |
// |                                           |
// =============================================

Test(round_up_div, u64)
{
	cr_assert(eq(u64, u64_round_up_div(0, 7), 0));
	cr_assert(eq(u64, u64_round_up_div(14, 7), 2));
	cr_assert(eq(u64, u64_round_up_div(15, 7), 3));
	cr_assert(eq(u64, u64_round_up_div(1, 1), 1));
	cr_assert(eq(u64, u64_round_up_div(63, 64), 1));
	// No overflow close to the maximum
	cr_assert(eq(u64, u64_round_up_div(UINT64_MAX, 2), UINT64_C(1) << 63));
	cr_assert(eq(u64, u64_round_up_div(UINT64_MAX, UINT64_MAX), 1));
	cr_assert(eq(u64, u64_round_up_div(UINT64_MAX - 1, UINT64_MAX), 1));
}

Test(round_up_div, i32_all_signs)
{
	cr_assert(eq(i32, i32_round_up_div(0, 3), 0));
	cr_assert(eq(i32, i32_round_up_div(7, 3), 3));
	cr_assert(eq(i32, i32_round_up_div(6, 3), 2));
	// ceil of negative quotients
	cr_assert(eq(i32, i32_round_up_div(-7, 3), -2));
	cr_assert(eq(i32, i32_round_up_div(-6, 3), -2));
	cr_assert(eq(i32, i32_round_up_div(-4, 2), -2));
	cr_assert(eq(i32, i32_round_up_div(7, -3), -2));
	cr_assert(eq(i32, i32_round_up_div(-7, -3), 3));
	// No overflow close to the maximum
	cr_assert(eq(i32, i32_round_up_div(INT32_MAX, 2), (INT32_MAX / 2) + 1));
	cr_assert(eq(i32, i32_round_up_div(INT32_MAX, INT32_MAX), 1));
}

Test(round_up_div, generic_macro_dispatch)
{
	uint64_t u = 10;
	int32_t i  = -10;
	cr_assert(eq(u64, INT_ROUND_UP_DIV(u, (uint64_t)4), 3));
	cr_assert(eq(i32, INT_ROUND_UP_DIV(i, 4), -2));
}

// =============================================
// |                                           |
// |               next_pow2_log               |
// |                                           |
// =============================================

Test(next_pow2_log, values)
{
	cr_assert(eq(u64, next_pow2_log(0), 0));
	cr_assert(eq(u64, next_pow2_log(1), 0));
	cr_assert(eq(u64, next_pow2_log(2), 1));
	cr_assert(eq(u64, next_pow2_log(3), 2));
	cr_assert(eq(u64, next_pow2_log(4), 2));
	cr_assert(eq(u64, next_pow2_log(5), 3));
	for (int b = 1; b < 64; ++b)
	{
		cr_assert(eq(u64, next_pow2_log(UINT64_C(1) << b), (uint64_t)b));
		cr_assert(eq(u64, next_pow2_log((UINT64_C(1) << b) + 1), (uint64_t)b + 1));
	}
	cr_assert(eq(u64, next_pow2_log(UINT64_MAX), 64));
}
