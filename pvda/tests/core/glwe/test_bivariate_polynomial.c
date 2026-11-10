#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/parameterized.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <criterion/redirect.h>

#include "backend_arithmetic.h"

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "logger.h"
#include "math_utils.h"
#include "random_polynomial.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

//! COMMON PART (begin)

PvdaParamTest(poly_biv_size, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	cr_assert(eq(i64, glwe_params_n_limbs(params_glwe), params_glwe->ciphertext_nb_limbs));

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(biv_bytes, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	cr_assert(eq(i64, biv_bytes(params_glwe), params_glwe->nn * glwe_params_l_a(params_glwe) * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(biv_to_univ_rnx, runs, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyBiv* pol                   = new_biv(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	//Set the first limb for the first coefficient to 1, ie, first coefficient is set to 2^-kappa
	pol->ptr[0] = 1;

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol, 0);

	//Checks that the first coefficient is indeed 2^-kappa
	cr_assert(eq(dbl, pol_univ_computed[0], ldexp(1.0, -params_glwe->kappa)));

	//And the rest should be 0
	for (uint64_t p = 1; p < params_glwe->nn; p++)
		cr_assert(eq(dbl, pol_univ_computed[p], 0), "pol_univ_computed[%ld] = %ld ", p, pol_univ_computed[p]);

	delete_biv(pol);
	delete_univ_rnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_rnx_to_biv, one_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv(params_glwe);

	memset(pol_univ, 0, univ_rnx_bytes(params_glwe));
	// Define pol_univ(X) = 2^{-params->kappa} , i.e. in Zn[X,Y] pol(X,Y) = Y
	pol_univ[0] = ldexp(1.0, -params_glwe->kappa);

	// Compute pol_univ's base-2params->kappa normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	// Asserts pol_computed = Y
	cr_assert(eq(int, pol_computed->ptr[0], 1), "pol_computed[%ld, %ld] = %ld ", 0, 1, pol_computed[0]);
	for (uint64_t i = 1; i < glwe_params_l_a(params_glwe) * params_glwe->nn; i++)
		cr_assert(eq(int, pol_computed->ptr[i], 0), "pol_computed[%ld, %ld] = %ld ", i / params_glwe->nn,
		          i % params_glwe->nn, pol_computed[i]);

	delete_univ_rnx(pol_univ);
	delete_biv(pol_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_rnx_to_biv, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv(params_glwe);

	// Draws normaly pol_univ in Rn[X] (!= torus)
	rnx_random_vec(module, pol_univ, params_glwe);

	// Computes pol_univ's base-2params->kappa normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	double err_length = glwe_bivariate_epsilon(params_glwe) + 3 * DBL_EPSILON;

	// Asserts pol_computed, in Rn[X] (with Y = 2^{-params->kappa}), is equal to pol_univ
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double pol_computed_p = 0;
		for (uint64_t i = 1; i <= glwe_params_l_a(params_glwe); i++)
			pol_computed_p += ldexp((double)pol_computed->ptr[(i - 1) * params_glwe->nn + p], -i * params_glwe->kappa);

		cr_assert(epsilon_eq(dbl, rnx_torus_distance(pol_computed_p, pol_univ[p]), 0, err_length));
	}

	delete_univ_rnx(pol_univ);
	delete_biv(pol_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

// univ_rnx_to_biv must overwrite every limb of the output, whatever the buffer held before.
// Converts the same input into a zeroed buffer and into a garbage-filled buffer: both results must be identical.
// The input sweeps powers of two from 2^-1 down to 2^-(kappa * l_a), so that small values (which have fewer
// significant limbs than the depth allows) are covered.
PvdaParamTest(univ_rnx_to_biv, overwrites_reused_buffer, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	uint64_t nn       = params_glwe->nn;
	uint64_t l_a      = glwe_params_l_a(params_glwe);
	uint64_t max_bits = params_glwe->kappa * l_a;

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_zeroed   = new_biv(params_glwe);
	PolyBiv* pol_garbage  = new_biv(params_glwe);

	for (uint64_t p = 0; p < nn; p++)
	{
		double sign = (p & 1) ? -1.0 : 1.0;
		pol_univ[p] = sign * ldexp(1.0, -(int)(1 + p % max_bits));
	}

	for (uint64_t i = 0; i < l_a * nn; i++) pol_garbage->ptr[i] = 12345;

	univ_rnx_to_biv(params_glwe, pol_zeroed, pol_univ, 0);
	univ_rnx_to_biv(params_glwe, pol_garbage, pol_univ, 0);

	for (uint64_t i = 0; i < l_a * nn; i++)
		cr_assert(eq(i64, pol_garbage->ptr[i], pol_zeroed->ptr[i]), "limb %" PRIu64 ", coefficient %" PRIu64,
		          i / nn, i % nn);

	delete_univ_rnx(pol_univ);
	delete_biv(pol_zeroed);
	delete_biv(pol_garbage);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_rnx_to_biv, maths_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ          = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed          = new_biv(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	rnx_random_vec(module, pol_univ, params_glwe);

	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol_computed, 0);

	double err_length = glwe_params_l_a(params_glwe) + 3 * DBL_EPSILON;

	// Asserts pol_univ_computed(X) = pol_univ(X)
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		cr_assert(epsilon_eq(dbl, rnx_torus_distance(pol_univ[p], pol_univ_computed[p]), 0, err_length));
	}

	delete_univ_rnx(pol_univ);
	delete_biv(pol_computed);
	delete_univ_rnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_tnx_to_biv, maths_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	univ_znx_random_uniform_pow2(module, (int64_t*)pol_univ, params_glwe->nn, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed, 0);

	int bits = glwe_params_l_a(params_glwe) * params_glwe->kappa;
	assert_tnx_close_enough_vec(pol_univ, pol_univ_computed, params_glwe->nn, bits);

	delete_univ_tnx(pol_univ);
	delete_biv(pol_computed);
	delete_univ_tnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_tnx_rnx_to_biv, maths_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ_rnx      = new_univ_rnx(params_glwe);
	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	univ_znx_random_uniform_pow2(module, (int64_t*)pol_univ, params_glwe->nn, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed, 0);

	int bits = glwe_params_l_a(params_glwe) * params_glwe->kappa;
	assert_tnx_close_enough_vec(pol_univ, pol_univ_computed, params_glwe->nn, bits);

	delete_univ_tnx(pol_univ);
	delete_biv(pol_computed);
	delete_univ_tnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_tnx_to_biv, small_znx, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	univ_znx_random_uniform_pow2(module, (int64_t*)pol_univ, params_glwe->nn, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed, 0);

	int bits = glwe_params_l_a(params_glwe) * params_glwe->kappa;
	assert_tnx_close_enough_vec(pol_univ, pol_univ_computed, params_glwe->nn, bits);

	delete_univ_tnx(pol_univ);
	delete_biv(pol_computed);
	delete_univ_tnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(tnx_rnx_encoding, back_and_forth_tnx_via_biv, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	int vec_size = params_glwe->nn;

	PolyUnivRnX* rnx_values = new_univ_rnx(params_glwe);
	PolyUnivTnX* tnx_values = new_univ_tnx(params_glwe);
	PolyUnivTnX* tnx_final  = new_univ_tnx(params_glwe);
	PolyBiv* biv            = new_biv(params_glwe);

	int l     = glwe_params_l_a(params_glwe);
	int kappa = params_glwe->kappa;
	// Max precision according to biv representation
	uint64_t precisionk = l * kappa;

	// Max precision is limited by float precision
	uint64_t precision = precisionk > 53 ? 53 : precisionk;

	univ_znx_random_uniform_pow2(module, tnx_values, vec_size, 64);

	univ_tnx_to_biv(params_glwe, biv, tnx_values, 0);
	biv_to_univ_rnx(params_glwe, rnx_values, biv, 0);
	univ_rnx_to_tnx(params_glwe, tnx_final, rnx_values);

	assert_tnx_close_enough_vec(tnx_final, tnx_values, vec_size, precision);

	delete_univ_rnx(rnx_values);
	delete_univ_tnx(tnx_values);
	delete_univ_tnx(tnx_final);
	delete_biv(biv);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(tnx_rnx_encoding, back_and_forth_tnx_via_biv_minus_offset, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	int vec_size = params_glwe->nn;

	PolyUnivRnX* rnx_values = new_univ_rnx(params_glwe);
	PolyUnivTnX* tnx_values = new_univ_tnx(params_glwe);
	PolyUnivTnX* tnx_final  = new_univ_tnx(params_glwe);
	PolyBiv* biv            = new_biv(params_glwe);

	for (int i = -5; i < 5; ++i)
	{
		int l     = glwe_params_l_a(params_glwe);
		int kappa = params_glwe->kappa;

		// Max precision according to biv representation
		uint64_t precisionk = l * kappa;
		if (i > 0) precisionk -= i;

		// Max precision is limited by float precision
		uint64_t precision = precisionk >= 53 ? 53 : precisionk;

		univ_znx_random_uniform_pow2(module, tnx_values, vec_size, 56);

		univ_tnx_to_biv(params_glwe, biv, tnx_values, i);
		biv_to_univ_rnx(params_glwe, rnx_values, biv, i);
		univ_rnx_to_tnx(params_glwe, tnx_final, rnx_values);

		assert_tnx_close_enough_vec(tnx_final, tnx_values, vec_size, precision);
	}

	delete_univ_rnx(rnx_values);
	delete_univ_tnx(tnx_values);
	delete_univ_tnx(tnx_final);
	delete_biv(biv);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(tnx_rnx_encoding, back_and_forth_rnx_via_biv, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);
	int vec_size = params_glwe->nn;

	double biv_err = glwe_bivariate_epsilon(params_glwe);

	PolyUnivRnX* rnx_values = new_univ_rnx(params_glwe);
	PolyUnivTnX* tnx_values = new_univ_tnx(params_glwe);
	PolyUnivRnX* rnx_final  = new_univ_rnx(params_glwe);
	PolyBiv* biv            = new_biv(params_glwe);

	rnx_random_vec(module, rnx_values, params_glwe);
	for (int i = 0; i < vec_size; ++i)
	{
		rnx_values[i] = rnx_values[i] - floor(rnx_values[i]);
	}

	for (int i = 0; i < 5; ++i)
	{
		univ_rnx_to_biv(params_glwe, biv, rnx_values, i);
		biv_to_univ_tnx(params_glwe, tnx_values, biv, i);
		univ_tnx_to_rnx(params_glwe, rnx_final, tnx_values);

		double adj_biv_err = biv_err * (1 << i);

		pvda_assert_polynomial_distance(params_glwe, rnx_final, rnx_values, 3 * DBL_EPSILON + adj_biv_err,
		                                5 * DBL_EPSILON + adj_biv_err);
	}

	delete_biv(biv);
	delete_univ_rnx(rnx_values);
	delete_univ_tnx(tnx_values);
	delete_univ_rnx(rnx_final);

	DELETE_PVDA_PARAMS_GLWE;
}

//! BIV POLY PART (begin)

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
PvdaParamTest(biv_coef_number, classic_params, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Asserts biv_coef_number returns params_glwe->nn * params_glwe->l
	cr_assert(eq(i64, biv_coef_number(params_glwe), params_glwe->nn * glwe_params_l_a(params_glwe)));

	DELETE_PVDA_PARAMS_GLWE;
}

/*
// TEST DISABLED UNTIL NORMAL RNG IS REIMPLEMENTED

PvdaParamTest(biv_random_normal, does_not_crash, default_params_fn)
{
    INIT_PVDA_PARAMS_GLWE(param);

    PolyBiv* pol = new_biv(params_glwe);

    int status = biv_random_normal(params_glwe, pol);

    cr_assert(eq(int, status, 0, "normal_random_biv failed."));

    delete_biv(pol);

    DELETE_PVDA_PARAMS_GLWE;
}*/

/*
// TEST DISABLED UNTIL NORMAL RNG IS REIMPLEMENTED
PvdaParamTest(biv_random_normal, output_is_normalized, default_params_fn)
{
    INIT_PVDA_PARAMS_GLWE(param);

    PolyBiv* pol = new_biv(params_glwe);

    // Draw normaly pol in Zn[X,Y]
    biv_random_normal(params_glwe, pol);

    // Asserts pol is normalized.
    // i.e. that each coefficient is between -2^(params->kappa-1) (inclusive) and 2^(params->kappa-1) (exculsive)
    for (uint64_t i = 0; i < glwe_params_l_a(params_glwe) * params_glwe->nn; i++)
    {
        cr_assert(lt(i64, pol->ptr[i], (1LL << (params_glwe->kappa - 1))));
        cr_assert(ge(i64, pol->ptr[i], -(1LL << (params_glwe->kappa - 1))));
    }

    delete_biv(pol);

    DELETE_PVDA_PARAMS_GLWE;
}*/

/**
 * @brief Test biv_add correctness with random normal polynomials
 */
PvdaParamTest(biv_add, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyBiv* pol_lhs      = new_biv(params_glwe);
	PolyBiv* pol_rhs      = new_biv(params_glwe);
	PolyBiv* sum_observed = new_biv(params_glwe);

	biv_random_uniform(module, params_glwe, pol_lhs, glwe_params_l_a(params_glwe));
	biv_random_uniform(module, params_glwe, pol_rhs, glwe_params_l_a(params_glwe));

	biv_add(module, sum_observed, pol_lhs, pol_rhs);

	// Asserts sum_computed = pol_lhs + pol_rhs
	for (uint64_t i = 0; i < glwe_params_l_a(params_glwe) * params_glwe->nn; i++)
	{
		cr_assert(eq(i64, sum_observed->ptr[i], pol_lhs->ptr[i] + pol_rhs->ptr[i]));
	}

	delete_biv(pol_lhs);
	delete_biv(pol_rhs);
	delete_biv(sum_observed);

	DELETE_PVDA_PARAMS_GLWE;
}

// =============================================
// |                                           |
// |     Small fixed parameters (N = 64)       |
// |                                           |
// =============================================

// K = 12, l_a = 4 (48 bits of precision), l_b = 3
#define SN 64
#define SK 12

static GLWEParams* small_params(void) { return new_glwe_params(SN, 1, SK, 7, 0, NOISE_UNIFORM_POWER_OF_TWO); }

/** @brief A backend with a seeded RNG, so that the tests below are reproducible. */
static PvdaBackend* small_module(void)
{
	PvdaBackend* module = pvda_new_spqlios_backend(SN);
	PvdaSeed seed;
	for (int i = 0; i < PVDA_SEED_BYTES; ++i) seed.bytes[i] = (uint8_t)(i * 7 + 2);
	cr_assert(pvda_backend_set_rng(module, PVDA_RNG_CPU, &seed) == 0);
	return module;
}

static void assert_digits_in(const PolyBiv* pol, int64_t lo, int64_t hi)
{
	for (uint64_t l = 0; l < pol->nb_limbs; ++l)
		for (uint64_t p = 0; p < pol->nn; ++p)
		{
			int64_t d = pol->ptr[l * pol->stride + p];
			cr_assert(d >= lo && d <= hi, "limb %lu coef %lu: %ld not in [%ld, %ld]", l, p, d, lo, hi);
		}
}

// =============================================
// |                                           |
// |           Allocation and views            |
// |                                           |
// =============================================

Test(biv_alloc, new_biv_and_custom_params)
{
	GLWEParams* params_glwe = small_params();
	PolyBiv* a              = new_biv(params_glwe);
	cr_assert_not_null(a);
	cr_assert(eq(u64, a->nn, SN));
	cr_assert(eq(u64, a->nb_limbs, glwe_params_l_a(params_glwe)));
	cr_assert(eq(i64, a->stride, SN));
	for (uint64_t i = 0; i < biv_coef_number(params_glwe); ++i) cr_assert(eq(i64, a->ptr[i], 0));

	PolyBiv* b = new_biv_custom_params(32, 9);
	cr_assert_not_null(b);
	cr_assert(eq(u64, b->nn, 32));
	cr_assert(eq(u64, b->nb_limbs, 9));
	cr_assert(eq(i64, b->stride, 32));
	for (uint64_t i = 0; i < 32 * 9; ++i) cr_assert(eq(i64, b->ptr[i], 0));

	delete_biv(a);
	delete_biv(b);
	delete_biv(NULL);
	delete_glwe_params(params_glwe);
}

Test(biv_alloc, view_does_not_copy)
{
	int64_t buf[3 * 10] = {0};
	PolyBiv v           = biv_view(8, 3, 10, buf);
	cr_assert(eq(u64, v.nn, 8));
	cr_assert(eq(u64, v.nb_limbs, 3));
	cr_assert(eq(i64, v.stride, 10));
	cr_assert(eq(ptr, v.ptr, buf));
	v.ptr[2 * v.stride + 7] = 5;
	cr_assert(eq(i64, buf[27], 5));
}

Test(biv_sizes, coef_numbers_and_bytes)
{
	GLWEParams* params_glwe = small_params();
	uint64_t l_a            = glwe_params_l_a(params_glwe);
	cr_assert(eq(u64, biv_coef_number(params_glwe), l_a * SN));
	cr_assert(eq(u64, biv_dft_coef_number(params_glwe), l_a * SN / 2));
	cr_assert(eq(u64, biv_bytes(params_glwe), l_a * SN * sizeof(int64_t)));
	delete_glwe_params(params_glwe);
}

// =============================================
// |                                           |
// |                Operations                 |
// |                                           |
// =============================================

Test(biv_ops, add_sub_including_missing_limbs_and_strides)
{
	PvdaBackend* module = small_module();
	int64_t a[4 * SN], b[2 * SN], res[4 * 2 * SN];
	cr_assert(univ_znx_random_uniform_pow2(module, a, 4 * SN, 40) == 0);
	cr_assert(univ_znx_random_uniform_pow2(module, b, 2 * SN, 40) == 0);
	PolyBiv va = biv_view(SN, 4, SN, a);
	PolyBiv vb = biv_view(SN, 2, SN, b);
	// Strided result: the other half of every row must not be touched
	for (int i = 0; i < 4 * 2 * SN; ++i) res[i] = 99;
	PolyBiv vr = biv_view(SN, 4, 2 * SN, res);

	cr_assert(eq(int, biv_add(module, &vr, &va, &vb), 0));
	for (int l = 0; l < 4; ++l)
		for (int p = 0; p < SN; ++p)
		{
			int64_t expected = a[l * SN + p] + (l < 2 ? b[l * SN + p] : 0);
			cr_assert(eq(i64, res[l * 2 * SN + p], expected));
			cr_assert(eq(i64, res[l * 2 * SN + SN + p], 99));
		}

	cr_assert(eq(int, biv_sub(module, &vr, &va, &vb), 0));
	for (int l = 0; l < 4; ++l)
		for (int p = 0; p < SN; ++p)
			cr_assert(eq(i64, res[l * 2 * SN + p], a[l * SN + p] - (l < 2 ? b[l * SN + p] : 0)));

	// b - a: limbs missing in b are read as 0
	cr_assert(eq(int, biv_sub(module, &vr, &vb, &va), 0));
	for (int l = 0; l < 4; ++l)
		for (int p = 0; p < SN; ++p)
			cr_assert(eq(i64, res[l * 2 * SN + p], (l < 2 ? b[l * SN + p] : 0) - a[l * SN + p]));

	// In place: a - a = 0
	cr_assert(eq(int, biv_sub(module, &va, &va, &va), 0));
	for (int i = 0; i < 4 * SN; ++i) cr_assert(eq(i64, a[i], 0));

	pvda_delete_backend(module);
}

/** @brief res = a(X^p) mod X^n + 1, naively. */
static void naive_automorphism(int64_t* res, const int64_t* a, int64_t p, uint64_t n)
{
	memset(res, 0, n * sizeof(int64_t));
	uint64_t two_n = 2 * n;
	for (uint64_t i = 0; i < n; ++i)
	{
		uint64_t e = (uint64_t)((((__int128)i * p) % (__int128)two_n + two_n) % two_n);
		if (e < n) res[e] += a[i];
		else res[e - n] -= a[i];
	}
}

Test(biv_ops, automorphism_matches_naive)
{
	PvdaBackend* module = small_module();
	int64_t a[3 * SN], res[3 * SN], expected[SN];
	PolyBiv va = biv_view(SN, 3, SN, a);
	PolyBiv vr = biv_view(SN, 3, SN, res);
	cr_assert(univ_znx_random_uniform_pow2(module, a, 3 * SN, 20) == 0);

	static const int64_t ps[] = {1, 3, 5, -1, 2 * SN - 1, 2 * SN + 3};
	for (size_t k = 0; k < sizeof(ps) / sizeof(ps[0]); ++k)
	{
		cr_assert(eq(int, biv_automorphism(module, ps[k], &vr, &va), 0), "p = %ld", ps[k]);
		for (int l = 0; l < 3; ++l)
		{
			naive_automorphism(expected, a + l * SN, ps[k], SN);
			for (int p = 0; p < SN; ++p) cr_assert(eq(i64, res[l * SN + p], expected[p]), "p = %ld limb %d coef %d", ps[k], l, p);
		}
	}
	pvda_delete_backend(module);
}

Test(biv_ops, normalize_keeps_the_value_and_balances_digits)
{
	PvdaBackend* module     = small_module();
	GLWEParams* params_glwe = small_params();
	PolyBiv* a              = new_biv(params_glwe);
	PolyBiv* n              = new_biv(params_glwe);
	PolyUnivTnX before[SN], after[SN];

	// Digits far outside [-2^(K-1), 2^(K-1)): carries must propagate across limbs (and wrap mod 1)
	cr_assert(biv_random_uniform_pow2(module, a, 30) == 0);
	a->ptr[0] = INT64_C(1) << 40;  // pure integer part on limb 0: vanishes mod 1
	biv_to_univ_tnx(params_glwe, before, a, 0);

	cr_assert(eq(int, biv_normalize(module, SK, n, a), 0));
	assert_digits_in(n, -(INT64_C(1) << (SK - 1)), (INT64_C(1) << (SK - 1)) - 1);
	biv_to_univ_tnx(params_glwe, after, n, 0);
	for (int p = 0; p < SN; ++p) cr_assert(eq(u64, after[p], before[p]), "coef %d", p);

	// Normalizing a normalized polynomial is the identity, in place too
	cr_assert(eq(int, biv_normalize(module, SK, n, n), 0));
	biv_to_univ_tnx(params_glwe, after, n, 0);
	for (int p = 0; p < SN; ++p) cr_assert(eq(u64, after[p], before[p]));

	delete_biv(a);
	delete_biv(n);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |          Univ -> Biv edge cases           |
// |                                           |
// =============================================

Test(univ_znx_to_biv, value_is_m_over_2_pow_offset)
{
	GLWEParams* params_glwe = small_params();
	PolyBiv* res            = new_biv(params_glwe);
	PolyUnivZnx m[SN];
	PolyUnivTnX t[SN];
	for (int p = 0; p < SN; ++p) m[p] = (p % 2 ? -1 : 1) * ((INT64_C(0x123456789) * (p + 1)) % (INT64_C(1) << 40));
	m[0] = 0;
	m[1] = -1;
	m[2] = (INT64_C(1) << 40) - 1;

	// Exact as long as the offset fits in the K * l_a = 48 bits of precision
	static const int64_t offsets[] = {1, 12, 13, 47, 48};
	for (size_t k = 0; k < sizeof(offsets) / sizeof(offsets[0]); ++k)
	{
		cr_assert(eq(int, univ_znx_to_biv(params_glwe, res, m, offsets[k]), 0));
		assert_digits_in(res, -(INT64_C(1) << (SK - 1)), INT64_C(1) << (SK - 1));
		biv_to_univ_tnx(params_glwe, t, res, 0);
		for (int p = 0; p < SN; ++p)
			cr_assert(eq(u64, t[p], (uint64_t)m[p] << (64 - offsets[k])), "offset %ld coef %d", offsets[k], p);
	}

	// A ZnX polynomial is 0 on the torus: non-positive offsets give 0 (and overwrite previous content)
	set_log_flag(LOG_ERROR);
	for (int64_t off = 0; off >= -3; --off)
	{
		res->ptr[5] = 77;
		cr_assert(eq(int, univ_znx_to_biv(params_glwe, res, m, off), 0));
		for (uint64_t i = 0; i < biv_coef_number(params_glwe); ++i) cr_assert(eq(i64, res->ptr[i], 0), "offset %ld", off);
	}

	delete_biv(res);
	delete_glwe_params(params_glwe);
}

Test(univ_rnx_to_biv, edge_values)
{
	GLWEParams* params_glwe = small_params();
	PolyBiv* res            = new_biv(params_glwe);
	PolyUnivRnX r[SN] = {0.0, -0.0, 3.0, -2.0, -0.25, 0.5, -0.5, 0x1p-48, -0x1p-48, 0x1p-60, 0.75, 1.25};
	// Expected torus values (truncated to 48 bits)
	const PolyUnivTnX expected[12] = {0, 0, 0, 0, UINT64_C(3) << 62, UINT64_C(1) << 63, UINT64_C(1) << 63,
	                                  UINT64_C(1) << 16, (uint64_t)-(INT64_C(1) << 16), 0, UINT64_C(3) << 62,
	                                  UINT64_C(1) << 62};
	PolyUnivTnX t[SN];

	for (uint64_t i = 0; i < biv_coef_number(params_glwe); ++i) res->ptr[i] = 1234;  // must be overwritten
	cr_assert(eq(int, univ_rnx_to_biv(params_glwe, res, r, 0), 0));
	assert_digits_in(res, -(INT64_C(1) << (SK - 1)), INT64_C(1) << (SK - 1));
	biv_to_univ_tnx(params_glwe, t, res, 0);
	for (int p = 0; p < 12; ++p) cr_assert(eq(u64, t[p], expected[p]), "coef %d (%a)", p, r[p]);
	for (int p = 12; p < SN; ++p) cr_assert(eq(u64, t[p], 0));

	delete_biv(res);
	delete_glwe_params(params_glwe);
}

Test(univ_tnx_to_biv, strided_output_and_offsets)
{
	GLWEParams* params_glwe = small_params();
	int64_t buf[4 * 2 * SN];
	PolyBiv res = biv_view(SN, 4, 2 * SN, buf);
	PolyUnivTnX in[SN], back[SN];
	for (int p = 0; p < SN; ++p) in[p] = UINT64_C(0x9e3779b97f4a7c15) * (uint64_t)(p + 1);
	in[0] = 0;
	in[1] = UINT64_C(1) << 63;
	in[2] = UINT64_MAX;

	for (int i = 0; i < 4 * 2 * SN; ++i) buf[i] = 99;
	cr_assert(eq(int, univ_tnx_to_biv(params_glwe, &res, in, 0), 0));
	for (int l = 0; l < 4; ++l)
		for (int p = 0; p < SN; ++p) cr_assert(eq(i64, buf[l * 2 * SN + SN + p], 99), "gap overwritten");
	assert_digits_in(&res, -(INT64_C(1) << (SK - 1)), INT64_C(1) << (SK - 1));

	// Back to TnX: equal up to the 48 bits of precision (the input LSB is dropped, then truncated)
	biv_to_univ_tnx(params_glwe, back, &res, 0);
	for (int p = 0; p < SN; ++p) cr_assert(le(u64, tnx_torus_distance(back[p], in[p]), UINT64_C(1) << 16), "coef %d", p);

	// bit_offset divides by 2^offset, and biv_to_univ_tnx with the same offset multiplies back
	cr_assert(eq(int, univ_tnx_to_biv(params_glwe, &res, in, 8), 0));
	biv_to_univ_tnx(params_glwe, back, &res, 8);
	for (int p = 0; p < SN; ++p) cr_assert(le(u64, tnx_torus_distance(back[p], in[p]), UINT64_C(1) << 24), "coef %d", p);

	delete_glwe_params(params_glwe);
}

Test(biv_to_univ_tnx, digit_weights)
{
	GLWEParams* params_glwe = small_params();
	PolyBiv* pol            = new_biv(params_glwe);
	PolyUnivTnX t[SN];
	// Coefficient p has a single digit 1 on limb p % 4: value 2^(-K (limb + 1))
	for (int p = 0; p < SN; ++p) pol->ptr[(p % 4) * SN + p] = (p / 4) % 2 ? -1 : 1;
	biv_to_univ_tnx(params_glwe, t, pol, 0);
	for (int p = 0; p < SN; ++p)
	{
		uint64_t w = UINT64_C(1) << (64 - SK * (p % 4 + 1));
		cr_assert(eq(u64, t[p], (p / 4) % 2 ? -w : w), "coef %d", p);
	}
	// Positive offset multiplies by 2^offset, negative divides
	biv_to_univ_tnx(params_glwe, t, pol, 4);
	cr_assert(eq(u64, t[1], UINT64_C(1) << (64 - 2 * SK + 4)));
	biv_to_univ_tnx(params_glwe, t, pol, -4);
	cr_assert(eq(u64, t[1], UINT64_C(1) << (64 - 2 * SK - 4)));
	// A digit shifted beyond the 64 bits disappears
	biv_to_univ_tnx(params_glwe, t, pol, -64);
	cr_assert(eq(u64, t[0], 0));
	delete_biv(pol);
	delete_glwe_params(params_glwe);
}

// =============================================
// |                                           |
// |          DFT and prepared forms           |
// |                                           |
// =============================================

Test(biv_dft, round_trip)
{
	PvdaBackend* module     = small_module();
	GLWEParams* params_glwe = small_params();
	PolyBiv* a              = new_biv(params_glwe);
	PolyBiv* back           = new_biv(params_glwe);
	BE_PolyBivDFT* dft      = new_biv_dft(module, params_glwe);
	BE_PolyBivDFT* dft2     = new_biv_dft_custom_params(module, 7);
	cr_assert_not_null(dft);
	cr_assert_not_null(dft2);

	cr_assert(biv_random_uniform_pow2(module, a, 30) == 0);
	cr_assert(eq(int, biv_to_dft(module, params_glwe, dft, a), 0));
	cr_assert(eq(int, biv_from_dft(module, params_glwe, back, dft), 0));
	for (uint64_t i = 0; i < biv_coef_number(params_glwe); ++i) cr_assert(eq(i64, back->ptr[i], a->ptr[i]));

	delete_biv_dft(module, dft2);
	delete_biv_dft(module, dft);
	delete_biv(back);
	delete_biv(a);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

Test(biv_prep, prepared_vector_gives_the_same_vmp)
{
	PvdaBackend* module     = small_module();
	GLWEParams* params_glwe = small_params();
	const uint64_t nrows = glwe_params_l_a(params_glwe), ncols = 2;

	PolyBiv* a               = new_biv(params_glwe);
	PolyBiv* mat             = new_biv_custom_params(SN, nrows * ncols);
	BE_PolyMatVmpPrep* pmat  = pvda_new_vmp_pmat(module, nrows, ncols);
	BE_PolyBivVmpPrep* a_prep = new_biv_prep(module, params_glwe);
	BE_PolyBivVmpPrep* other  = new_biv_prep_custom_params(module, 3);
	BE_PolyBivDFT* r1        = pvda_new_biv_dft(module, ncols);
	BE_PolyBivDFT* r2        = pvda_new_biv_dft(module, ncols);
	PolyBiv* c1              = new_biv_custom_params(SN, ncols);
	PolyBiv* c2              = new_biv_custom_params(SN, ncols);
	cr_assert_not_null(a_prep);
	cr_assert_not_null(other);

	cr_assert(biv_random_uniform_pow2(module, a, 12) == 0);
	cr_assert(biv_random_uniform_pow2(module, mat, 12) == 0);
	cr_assert(eq(int, pvda_biv_to_vmp_pmat(module, pmat, mat, nrows, ncols), 0));

	cr_assert(eq(int, pvda_vmp_apply_dft(module, r1, ncols, a, pmat, nrows, ncols), 0));
	cr_assert(eq(int, biv_to_prep(module, params_glwe, a_prep, a), 0));
	cr_assert(eq(int, pvda_vmp_apply_prepared_to_dft(module, r2, ncols, a_prep, nrows, pmat, nrows, ncols), 0));
	cr_assert(eq(int, pvda_biv_idft(module, c1, r1, ncols), 0));
	cr_assert(eq(int, pvda_biv_idft(module, c2, r2, ncols), 0));
	for (uint64_t i = 0; i < ncols * SN; ++i) cr_assert(eq(i64, c1->ptr[i], c2->ptr[i]), "coef %lu", i);

	delete_biv(c2);
	delete_biv(c1);
	pvda_delete_biv_dft(module, r2);
	pvda_delete_biv_dft(module, r1);
	delete_biv_prep(module, other);
	delete_biv_prep(module, a_prep);
	pvda_delete_vmp_pmat(module, pmat);
	delete_biv(mat);
	delete_biv(a);
	delete_glwe_params(params_glwe);
	pvda_delete_backend(module);
}

// =============================================
// |                                           |
// |                  Prints                   |
// |                                           |
// =============================================

Test(biv_print, printf_poly_biv, .init = cr_redirect_stdout)
{
	int64_t buf[4] = {1, -2, 3, 4};
	PolyBiv v      = biv_view(2, 2, 2, buf);
	printf_poly_biv(&v);
	fflush(stdout);
	cr_assert_stdout_eq_str("\n\nY^0 : 1 X^0 -2 X^1 \nY^1 : 3 X^0 4 X^1 ");
}

Test(biv_print, print_coefs_biv_limits, .init = cr_redirect_stdout)
{
	int64_t buf[6] = {1, -2, 30, 4, 5, 6};
	PolyBiv v      = biv_view(2, 2, 3, buf);  // strided: buf[2] and buf[5] are not part of the polynomial
	print_coefs_biv(&v, 1, 5);                // at most 1 coefficient, limbs capped by nb_limbs
	print_coefs_biv(&v, 10, 1);               // coefficients capped by nn
	fflush(stdout);
	cr_assert_stdout_eq_str("Biv lvl 000: 001 \nBiv lvl 001: 004 \nBiv lvl 000: 001 -02 \n");
}
