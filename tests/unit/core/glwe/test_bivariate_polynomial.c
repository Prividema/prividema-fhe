#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/parameterized.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! COMMON PART (begin)

PvdaParamTest(poly_biv_size, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	cr_assert(eq(i64, glwe_params_n_limbs(params_glwe), params_glwe->ciphertext_nb_limbs));

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(poly_biv_bytes, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	cr_assert(eq(i64, poly_biv_bytes(params_glwe), params_glwe->nn * glwe_params_l_a(params_glwe) * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(biv_to_univ_rnx, runs, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyBiv* pol                   = new_biv_poly(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	//Set the first limb for the first coefficient to 1, ie, first coefficient is set to 2^-kappa
	pol[0] = 1;

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol);

	//Checks that the first coefficient is indeed 2^-kappa
	cr_assert(eq(dbl, pol_univ_computed[0], ldexp(1.0, -params_glwe->kappa)));

	//And the rest should be 0
	for (uint64_t p = 1; p < params_glwe->nn; p++)
		cr_assert(eq(dbl, pol_univ_computed[p], 0), "pol_univ_computed[%ld] = %ld ", p, pol_univ_computed[p]);

	free(pol);
	delete_univ_rnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_rnx_to_biv, one_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv_poly(params_glwe);

	memset(pol_univ, 0, poly_univ_rnx_bytes(params_glwe));
	// Define pol_univ(X) = 2^{-params->kappa} , i.e. in Zn[X,Y] pol(X,Y) = Y
	pol_univ[0] = ldexp(1.0, -params_glwe->kappa);

	// Compute pol_univ's base-2params->kappa normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	// Asserts pol_computed = Y
	cr_assert(eq(int, pol_computed[0], 1), "pol_computed[%ld, %ld] = %ld ", 0, 1, pol_computed[0]);
	for (uint64_t i = 1; i < glwe_params_l_a(params_glwe) * params_glwe->nn; i++)
		cr_assert(eq(int, pol_computed[i], 0), "pol_computed[%ld, %ld] = %ld ", i / params_glwe->nn,
		          i % params_glwe->nn, pol_computed[i]);

	delete_univ_rnx(pol_univ);
	free(pol_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_rnx_to_biv, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv_poly(params_glwe);

	// Draws normaly pol_univ in Rn[X] (!= torus)
	normal_random_vec(pol_univ, params_glwe->nn, 0.0, 1e-2);

	// Computes pol_univ's base-2params->kappa normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	double err_length = glwe_bivariate_epsilon(params_glwe) + 3 * DBL_EPSILON;

	// Asserts pol_computed, in Rn[X] (with Y = 2^{-params->kappa}), is equal to pol_univ
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double pol_computed_p = 0;
		for (uint64_t i = 1; i <= glwe_params_l_a(params_glwe); i++)
			pol_computed_p += ldexp((double)pol_computed[(i - 1) * params_glwe->nn + p], -i * params_glwe->kappa);

		cr_assert(epsilon_eq(dbl, torus_distance(pol_computed_p, pol_univ[p]), 0, err_length));
	}

	delete_univ_rnx(pol_univ);
	free(pol_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_rnx_to_biv, maths_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ          = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	normal_random_vec(pol_univ, params_glwe->nn, 0.0, 1e-2);

	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ, 0);

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol_computed);

	double err_length = glwe_params_l_a(params_glwe) + 3 * DBL_EPSILON;

	// Asserts pol_univ_computed(X) = pol_univ(X)
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		cr_assert(epsilon_eq(dbl, torus_distance(pol_univ[p], pol_univ_computed[p]), 0, err_length));
	}

	delete_univ_rnx(pol_univ);
	free(pol_computed);
	delete_univ_rnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_tnx_to_biv, maths_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	uniform_random_vec(params_glwe->nn, (PolyUniv*)pol_univ, 1, params_glwe->nn, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed);

	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		int bits = glwe_params_l_a(params_glwe) * params_glwe->kappa;
		if (bits >= 64)
			cr_assert(eq(u64, pol_univ[p], pol_univ_computed[p]));
		else
		{
			uint64_t max_diff = 1ULL << (64 - bits);
			uint64_t diff     = pol_univ[p] > pol_univ_computed[p] ? pol_univ[p] - pol_univ_computed[p]
			                                                       : pol_univ_computed[p] - pol_univ[p];
			cr_assert(lt(u64, diff, max_diff));
		}
	}

	delete_univ_tnx(pol_univ);
	free(pol_computed);
	delete_univ_tnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_tnx_rnx_to_biv, maths_test, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivRnX* pol_univ_rnx      = new_univ_rnx(params_glwe);
	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	uniform_random_vec(params_glwe->nn, (PolyUniv*)pol_univ, 1, params_glwe->nn, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed);

	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		int bits = glwe_params_l_a(params_glwe) * params_glwe->kappa;
		if (bits >= 64)
			cr_assert(eq(u64, pol_univ[p], pol_univ_computed[p]));
		else
		{
			uint64_t max_diff = 1ULL << (64 - bits);
			uint64_t diff     = pol_univ[p] > pol_univ_computed[p] ? pol_univ[p] - pol_univ_computed[p]
			                                                       : pol_univ_computed[p] - pol_univ[p];
			cr_assert(lt(u64, diff, max_diff));
		}
	}

	delete_univ_tnx(pol_univ);
	free(pol_computed);
	delete_univ_tnx(pol_univ_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(univ_tnx_to_biv, small_znx, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	uniform_random_vec(params_glwe->nn, (PolyUniv*)pol_univ, 1, params_glwe->nn, 12);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed);

	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		int bits = glwe_params_l_a(params_glwe) * params_glwe->kappa;
		if (bits >= 64)
			cr_assert(eq(u64, pol_univ[p], pol_univ_computed[p]));
		else
		{
			uint64_t max_diff = 1ULL << (64 - bits);
			uint64_t diff     = pol_univ[p] > pol_univ_computed[p] ? pol_univ[p] - pol_univ_computed[p]
			                                                       : pol_univ_computed[p] - pol_univ[p];
			cr_assert(lt(u64, diff, max_diff));
		}
	}

	delete_univ_tnx(pol_univ);
	free(pol_computed);
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
	PolyBiv* biv            = new_biv_poly(params_glwe);

	uint64_t precision1 = (1ULL << (64 - 53));

	int l               = glwe_params_l_a(params_glwe);
	int kappa           = params_glwe->kappa;
	uint64_t precision2 = l * kappa >= 64 ? 0 : (1ULL << (64 - l * kappa));

	uint64_t precision = precision1 > precision2 ? precision1 : precision2;

	uniform_random_pol_znx(tnx_values, vec_size, 64);

	univ_tnx_to_biv(params_glwe, biv, tnx_values);
	biv_to_univ_rnx(params_glwe, rnx_values, biv);
	univ_rnx_to_tnx(params_glwe, tnx_final, rnx_values);

	for (uint64_t i = 0; i < vec_size; ++i)
	{
		uint64_t diff = tnx_final[i] > tnx_values[i] ? tnx_final[i] - tnx_values[i] : tnx_values[i] - tnx_final[i];
		cr_assert(le(u64, diff, precision));
	}

	delete_univ_rnx(rnx_values);
	delete_univ_tnx(tnx_values);
	delete_univ_tnx(tnx_final);
	free(biv);

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
	PolyBiv* biv            = new_biv_poly(params_glwe);

	normal_random_vec(rnx_values, vec_size, 0, 0.2);
	for (int i = 0; i < vec_size; ++i)
	{
		rnx_values[i] = rnx_values[i] - floor(rnx_values[i]);
	}
	univ_rnx_to_biv(params_glwe, biv, rnx_values, 0);
	biv_to_univ_tnx(params_glwe, tnx_values, biv);
	univ_tnx_to_rnx(params_glwe, rnx_final, tnx_values);

	pvda_assert_polynomial_distance(params_glwe, rnx_final, rnx_values, 3 * DBL_EPSILON + biv_err,
	                                4 * DBL_EPSILON + biv_err);

	free(biv);
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
PvdaParamTest(poly_biv_coef_number, classic_params, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Asserts poly_biv_coef_number returns params_glwe->nn * params_glwe->l
	cr_assert(eq(i64, poly_biv_coef_number(params_glwe), params_glwe->nn * glwe_params_l_a(params_glwe)));

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(normal_random_biv_poly, does_not_crash, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyBiv* pol = new_biv_poly(params_glwe);

	int status = normal_random_biv_poly(params_glwe, pol);

	cr_assert(eq(int, status, 0, "normal_random_biv failed."));

	free(pol);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(normal_random_biv_poly, output_is_normalized, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyBiv* pol = new_biv_poly(params_glwe);

	// Draw normaly pol in Zn[X,Y]
	normal_random_biv_poly(params_glwe, pol);

	// Asserts pol is normalized.
	// i.e. that each coefficient is between -2^(params->kappa-1) (inclusive) and 2^(params->kappa-1) (exculsive)
	for (uint64_t i = 0; i < glwe_params_l_a(params_glwe) * params_glwe->nn; i++)
	{
		cr_assert(lt(i64, pol[i], (1LL << (params_glwe->kappa - 1))));
		cr_assert(ge(i64, pol[i], -(1LL << (params_glwe->kappa - 1))));
	}

	free(pol);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
PvdaParamTest(add_biv_poly, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	PolyBiv* pol_lhs      = new_biv_poly(params_glwe);
	PolyBiv* pol_rhs      = new_biv_poly(params_glwe);
	PolyBiv* sum_observed = new_biv_poly(params_glwe);

	normal_random_biv_poly(params_glwe, pol_lhs);
	normal_random_biv_poly(params_glwe, pol_rhs);

	add_biv_poly(module, params_glwe, sum_observed, pol_lhs, pol_rhs);

	// Asserts sum_computed = pol_lhs + pol_rhs
	for (uint64_t i = 0; i < glwe_params_l_a(params_glwe) * params_glwe->nn; i++)
	{
		cr_assert(eq(i64, sum_observed[i], pol_lhs[i] + pol_rhs[i]));
	}

	free(pol_lhs);
	free(pol_rhs);
	free(sum_observed);

	DELETE_PVDA_PARAMS_GLWE;
}
