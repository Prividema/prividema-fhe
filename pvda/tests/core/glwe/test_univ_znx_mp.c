#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/parameterized.h>
#include <stdint.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "math_utils.h"
#include "random_polynomial.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

// Tests of the multi-precision univariate polynomial type PolyUnivZnxMP and of its conversions.

static inline uint64_t words_for_bits(uint64_t bits) { return (bits + 63) / 64; }

/** @brief Fills the first l_a limbs of pol with normalized random digits in [-2^(K-1), 2^(K-1)). */
static void random_normalized_biv(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* pol)
{
	for (uint64_t i = 0; i < glwe_params_l_a(params_glwe); ++i)
		cr_assert(univ_znx_random_uniform_pow2(module, pol->ptr + i * pol->stride, pol->nn, params_glwe->kappa) == 0);
}

// =============================================
// |                                           |
// |           Allocation and ZnX              |
// |                                           |
// =============================================

Test(univ_znx_mp, new_is_zeroed)
{
	PolyUnivZnxMP* pol = new_univ_znx_mp(16, 3);
	cr_assert_not_null(pol);
	cr_assert(eq(u64, pol->nn, 16));
	cr_assert(eq(u64, pol->nb_words, 3));
	cr_assert(eq(u64, univ_znx_mp_bytes(pol), 16 * 3 * sizeof(uint64_t)));
	for (uint64_t i = 0; i < 16 * 3; ++i) cr_assert(eq(u64, pol->ptr[i], 0));
	delete_univ_znx_mp(pol);
}

Test(univ_znx_mp, new_rejects_empty)
{
	cr_assert_null(new_univ_znx_mp(0, 2));
	cr_assert_null(new_univ_znx_mp(8, 0));
	delete_univ_znx_mp(NULL);  // must not crash
}

Test(univ_znx_mp, znx_round_trip_sign_extends)
{
	const uint64_t nn  = 6;
	PolyUnivZnx in[6]  = {0, 1, -1, INT64_MAX, INT64_MIN, -123456789};
	PolyUnivZnx out[6] = {0};

	PolyUnivZnxMP* mp = new_univ_znx_mp(nn, 3);
	cr_assert(univ_znx_to_znx_mp(mp, in) == 0);

	for (uint64_t p = 0; p < nn; ++p)
	{
		uint64_t ext = in[p] < 0 ? UINT64_MAX : 0;
		cr_assert(eq(u64, mp->ptr[p * 3], (uint64_t)in[p]));
		cr_assert(eq(u64, mp->ptr[p * 3 + 1], ext));
		cr_assert(eq(u64, mp->ptr[p * 3 + 2], ext));
	}

	cr_assert(univ_znx_mp_to_znx(out, mp) == 0);
	for (uint64_t p = 0; p < nn; ++p) cr_assert(eq(i64, out[p], in[p]));

	delete_univ_znx_mp(mp);
}

Test(univ_znx_mp, to_znx_detects_overflow)
{
	PolyUnivZnx out[2];
	PolyUnivZnxMP* mp = new_univ_znx_mp(2, 2);

	// 2^64: does not fit in an int64_t, truncated to 0
	mp->ptr[0] = 0;
	mp->ptr[1] = 1;
	// 2^63 (positive, upper word 0): does not fit either
	mp->ptr[2] = 1ULL << 63;
	mp->ptr[3] = 0;

	cr_assert(univ_znx_mp_to_znx(out, mp) == -1);
	cr_assert(eq(i64, out[0], 0));
	cr_assert(eq(i64, out[1], INT64_MIN));

	delete_univ_znx_mp(mp);
}

// =============================================
// |                                           |
// |              MP <-> Biv                   |
// |                                           |
// =============================================

PvdaParamTest(univ_znx_mp, single_limb_value, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	const uint64_t nn    = params_glwe->nn;
	const uint64_t kappa = params_glwe->kappa;
	const uint64_t nw    = words_for_bits(kappa * glwe_params_l_a(params_glwe));

	// P(X, Y) = Y, i.e. 2^-K as a torus element
	PolyBiv* pol = new_biv(params_glwe);
	pol->ptr[0]  = 1;

	PolyUnivZnxMP* mp = new_univ_znx_mp(nn, nw);
	cr_assert(biv_to_univ_znx_mp(params_glwe, mp, pol, (int64_t)(64 * nw)) == 0);

	// Expect 2^(64 nw - K) in coefficient 0 and 0 elsewhere
	uint64_t bit = 64 * nw - kappa;
	for (uint64_t w = 0; w < nw; ++w)
		cr_assert(eq(u64, mp->ptr[w], w == bit / 64 ? 1ULL << (bit % 64) : 0), "word %lu", w);
	for (uint64_t i = nw; i < nn * nw; ++i) cr_assert(eq(u64, mp->ptr[i], 0));

	// And back
	PolyBiv* back = new_biv(params_glwe);
	cr_assert(univ_znx_mp_to_biv(params_glwe, back, mp, (int64_t)(64 * nw)) == 0);
	cr_assert(eq(i64, back->ptr[0], 1));
	for (uint64_t i = 1; i < glwe_params_l_a(params_glwe) * nn; ++i) cr_assert(eq(i64, back->ptr[i], 0));

	delete_biv(pol);
	delete_biv(back);
	delete_univ_znx_mp(mp);
	DELETE_PVDA_PARAMS_GLWE;
}

// Normalized biv -> MP -> biv is exact as soon as the MP coefficients hold K * l_a bits.
PvdaParamTest(univ_znx_mp, biv_round_trip_exact, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	const uint64_t nn = params_glwe->nn;
	const uint64_t l  = glwe_params_l_a(params_glwe);
	const uint64_t nw = words_for_bits(params_glwe->kappa * l);

	PolyBiv* pol      = new_biv(params_glwe);
	PolyBiv* back     = new_biv(params_glwe);
	PolyUnivZnxMP* mp = new_univ_znx_mp(nn, nw);
	random_normalized_biv(module, params_glwe, pol);

	cr_assert(biv_to_univ_znx_mp(params_glwe, mp, pol, (int64_t)(64 * nw)) == 0);
	cr_assert(univ_znx_mp_to_biv(params_glwe, back, mp, (int64_t)(64 * nw)) == 0);

	for (uint64_t i = 0; i < l * nn; ++i)
		cr_assert(eq(i64, back->ptr[i], pol->ptr[i]), "limb %lu coef %lu", i / nn, i % nn);

	delete_biv(pol);
	delete_biv(back);
	delete_univ_znx_mp(mp);
	DELETE_PVDA_PARAMS_GLWE;
}

// MP -> biv -> MP with bit_offset = K * l_a keeps exactly the K * l_a low bits (a / 2^(K l_a) mod 1).
PvdaParamTest(univ_znx_mp, mp_round_trip_mod_one, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	const uint64_t nn   = params_glwe->nn;
	const uint64_t bits = params_glwe->kappa * glwe_params_l_a(params_glwe);
	const uint64_t nw   = words_for_bits(bits) + 1;  // one more word than needed: holds an integer part

	PolyUnivZnxMP* mp   = new_univ_znx_mp(nn, nw);
	PolyUnivZnxMP* back = new_univ_znx_mp(nn, nw);
	PolyBiv* pol        = new_biv(params_glwe);
	cr_assert(univ_znx_random_uniform_pow2(module, (int64_t*)mp->ptr, nn * nw, 64) == 0);

	cr_assert(univ_znx_mp_to_biv(params_glwe, pol, mp, (int64_t)bits) == 0);
	cr_assert(biv_to_univ_znx_mp(params_glwe, back, pol, (int64_t)bits) == 0);

	for (uint64_t p = 0; p < nn; ++p)
		for (uint64_t w = 0; w < nw; ++w)
		{
			uint64_t lo   = 64 * w;
			uint64_t keep = bits <= lo ? 0 : (bits - lo >= 64 ? UINT64_MAX : (1ULL << (bits - lo)) - 1);
			cr_assert(eq(u64, back->ptr[p * nw + w] & keep, mp->ptr[p * nw + w] & keep), "coef %lu word %lu", p, w);
		}

	delete_univ_znx_mp(mp);
	delete_univ_znx_mp(back);
	delete_biv(pol);
	DELETE_PVDA_PARAMS_GLWE;
}

// With one word and K * l_a <= 64, the MP torus representation matches biv_to_univ_tnx.
PvdaParamTest(univ_znx_mp, matches_tnx_on_one_word, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	const uint64_t nn = params_glwe->nn;
	if (params_glwe->kappa * glwe_params_l_a(params_glwe) > 64)
	{
		DELETE_PVDA_PARAMS_GLWE;
		cr_skip_test("biv_to_univ_tnx is not exact when K * l_a > 64");
	}

	PolyBiv* pol      = new_biv(params_glwe);
	PolyUnivTnX* tnx  = new_univ_tnx(params_glwe);
	PolyUnivZnxMP* mp = new_univ_znx_mp(nn, 1);
	random_normalized_biv(module, params_glwe, pol);

	cr_assert(biv_to_univ_tnx(params_glwe, tnx, pol, 0) == 0);
	cr_assert(biv_to_univ_znx_mp(params_glwe, mp, pol, 64) == 0);
	for (uint64_t p = 0; p < nn; ++p) cr_assert(eq(u64, mp->ptr[p], tnx[p]), "coef %lu", p);

	delete_biv(pol);
	delete_univ_tnx(tnx);
	delete_univ_znx_mp(mp);
	DELETE_PVDA_PARAMS_GLWE;
}

// Rounding: the bits below the last limb are rounded to the nearest, ties up.
Test(univ_znx_mp, to_biv_rounds_to_nearest)
{
	// K = 8, l_a = 2: 16 bits of precision. Read 3-word values as torus elements (bit_offset = 192).
	GLWEParams* params_glwe = new_glwe_params(16, 1, 8, 2 * 2, 0, NOISE_UNIFORM_POWER_OF_TWO);
	cr_assert(eq(u64, glwe_params_l_a(params_glwe), 2));

	PolyUnivZnxMP* mp = new_univ_znx_mp(16, 3);
	PolyBiv* pol      = new_biv(params_glwe);

	// coef 0: 2^-16 + 2^-17 (tie) -> rounds up to 2 * 2^-16 -> limbs (0, 2)
	mp->ptr[0 * 3 + 2] = (1ULL << 48) + (1ULL << 47);
	// coef 1: 2^-16 + 2^-17 - 2^-192 (just below the tie) -> 2^-16 -> limbs (0, 1)
	mp->ptr[1 * 3 + 0] = UINT64_MAX;
	mp->ptr[1 * 3 + 1] = UINT64_MAX;
	mp->ptr[1 * 3 + 2] = (1ULL << 48) + (1ULL << 47) - 1;
	// coef 2: -2^-16 -> limbs (0, -1)
	mp->ptr[2 * 3 + 0] = 0;
	mp->ptr[2 * 3 + 1] = 0;
	mp->ptr[2 * 3 + 2] = (uint64_t)(-(int64_t)(1ULL << 48));
	// coef 3: 1/2 - 2^-17 (rounds to 1/2) -> carries up to -2^-1 -> limbs (-128, 0)
	mp->ptr[3 * 3 + 2] = (1ULL << 63) - (1ULL << 47);

	cr_assert(univ_znx_mp_to_biv(params_glwe, pol, mp, 192) == 0);
	const uint64_t nn = 16;
	cr_assert(eq(i64, pol->ptr[0], 0));
	cr_assert(eq(i64, pol->ptr[nn + 0], 2));
	cr_assert(eq(i64, pol->ptr[1], 0));
	cr_assert(eq(i64, pol->ptr[nn + 1], 1));
	cr_assert(eq(i64, pol->ptr[2], 0));
	cr_assert(eq(i64, pol->ptr[nn + 2], -1));
	cr_assert(eq(i64, pol->ptr[3], -128));
	cr_assert(eq(i64, pol->ptr[nn + 3], 0));

	delete_univ_znx_mp(mp);
	delete_biv(pol);
	delete_glwe_params(params_glwe);
}

Test(univ_znx_mp, rejects_dimension_mismatch)
{
	GLWEParams* params_glwe = new_glwe_params(16, 1, 8, 2 * 2, 0, NOISE_UNIFORM_POWER_OF_TWO);
	PolyUnivZnxMP* mp       = new_univ_znx_mp(8, 1);  // nn != 16
	PolyBiv* pol            = new_biv(params_glwe);

	cr_assert(univ_znx_mp_to_biv(params_glwe, pol, mp, 64) == -1);
	cr_assert(biv_to_univ_znx_mp(params_glwe, mp, pol, 64) == -1);

	delete_univ_znx_mp(mp);
	delete_biv(pol);
	delete_glwe_params(params_glwe);
}
