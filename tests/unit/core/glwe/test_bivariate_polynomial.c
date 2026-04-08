#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "common/rng.h"
#include "common/spqlios_alias.h"
#include "univariate_polynomial.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  19
#define NLIMBSBASE (KBASE + 1) * 4
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -7

//! COMMON PART (begin)

Test(poly_biv_size, basic)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_biv_size returns LBASE
	cr_assert(eq(i64, glwe_params_l(params_glwe), LBASE, "poly_biv_size failed: got %" PRId64 ", expected %" PRId64,
	             poly_biv_size(params_glwe), LBASE));

	delete_glwe_params(params_glwe);
}

Test(poly_biv_bytes, basic)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_biv_bytes returns NBASE * LBASE * sizeof(int64_t)
	cr_assert(eq(i64, poly_biv_bytes(params_glwe), NBASE * LBASE * sizeof(int64_t),
	             "poly_biv_bytes failed: got %" PRId64 ", expected %" PRId64, poly_biv_bytes(params_glwe),
	             NBASE * LBASE * sizeof(int64_t)));

	delete_glwe_params(params_glwe);
}

Test(biv_to_univ_rnx, runs)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyBiv* pol                   = new_biv_poly(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	//Set the first limb for the first coefficient to 1, ie, first coefficient is set to 2^-kappa
	pol[0] = 1;

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol);

	//Checks that the first coefficient is indeed 2^-kappa
	cr_assert(eq(dbl, pol_univ_computed[0], ldexp(1.0, -KAPPABASE)));

	//And the rest should be 0
	for (uint64_t p = 1; p < NBASE; p++)
		cr_assert(eq(dbl, pol_univ_computed[p], 0), "pol_univ_computed[%ld] = %ld ", p, pol_univ_computed[p]);

	free(pol);
	delete_univ_rnx(pol_univ_computed);
	delete_glwe_params(params_glwe);
}

Test(univ_to_biv, one_test)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv_poly(params_glwe);

	// Define pol_univ(X) = 2^{-KAPPABASE} , i.e. in Zn[X,Y] pol(X,Y) = Y
	pol_univ[0] = ldexp(1.0, -KAPPABASE);

	// Compute pol_univ's base-2KAPPABASE normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ);

	// Asserts pol_computed = Y
	cr_assert(eq(int, pol_computed[0], 1), "pol_computed[%ld, %ld] = %ld ", 0, 1, pol_computed[0]);

	for (uint64_t p = 1; p < NBASE; p++)
		cr_assert(eq(int, pol_computed[p], 0), "pol_computed[%ld, %ld] = %ld ", p, 1, pol_computed[p]);

	for (uint64_t i = 2; i <= LBASE; i++)
		for (uint64_t p = 0; p < NBASE; p++)
			cr_assert(eq(int, pol_computed[(i - 1) * NBASE + p], 0), "pol_computed[%ld, %ld] = %ld ", p, i,
			          pol_computed[(i - 1) * NBASE + p]);

	// Clean up
	delete_univ_rnx(pol_univ);
	free(pol_computed);
	delete_glwe_params(params_glwe);
}

Test(univ_to_biv, basic)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyUnivRnX* pol_univ = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed = new_biv_poly(params_glwe);

	// Draws normaly pol_univ in Rn[X] (!= torus)
	// TODO: use a torus element instead of rnx?
	normal_random_vec(pol_univ, NBASE, 0.0, 1e-2);

	// Computes pol_univ's base-2KAPPABASE normalized decomposition
	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ);

	double pol_computed_p = 0;

	// Asserts pol_computed, in Rn[X] (with Y = 2^{-KAPPABASE}), is equal to pol_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{
		pol_computed_p = 0;
		for (uint64_t i = 1; i <= LBASE; i++)
			pol_computed_p += ldexp((double)pol_computed[(i - 1) * NBASE + p], -i * KAPPABASE);

		//TODO: double check this and if correct add comment (and split into multiple lines) as to what it is doing
		cr_assert(epsilon_eq(dbl, pol_computed_p - floor(pol_computed_p) - pol_univ[p] + floor(pol_univ[p]), 0,
		                     ldexp(1.0, -(LBASE - 1) * KAPPABASE) + 3 * DBL_EPSILON));
	}

	delete_univ_rnx(pol_univ);
	free(pol_computed);
	delete_glwe_params(params_glwe);
}

Test(univ_rnx_to_biv, maths_test)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyUnivRnX* pol_univ          = new_univ_rnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivRnX* pol_univ_computed = new_univ_rnx(params_glwe);

	normal_random_vec(pol_univ, NBASE, 0.0, 1e-2);

	univ_rnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_rnx(params_glwe, pol_univ_computed, pol_computed);

	// Asserts pol_univ_computed(X) = pol_univ(X)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		cr_assert(epsilon_eq(dbl, pol_univ[p] - floor(pol_univ[p]) - pol_univ_computed[p] + floor(pol_univ_computed[p]),
		                     0, ldexp(1.0, -(LBASE - 1) * KAPPABASE) + 3 * DBL_EPSILON),
		          "pol_univ[%ld] = %lf and pol_univ_computed[%ld] = %lf", p, pol_univ[p], p, pol_univ_computed[p]);
	}

	// Clean up
	delete_univ_rnx(pol_univ);
	free(pol_computed);
	delete_univ_rnx(pol_univ_computed);
	delete_glwe_params(params_glwe);
}

Test(univ_tnx_to_biv, maths_test)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyUnivTnX* pol_univ          = new_univ_tnx(params_glwe);
	PolyBiv* pol_computed          = new_biv_poly(params_glwe);
	PolyUnivTnX* pol_univ_computed = new_univ_tnx(params_glwe);

	uniform_random_vec(NBASE, (PolyUniv*)pol_univ, 1, NBASE, 64);

	univ_tnx_to_biv(params_glwe, pol_computed, pol_univ);

	biv_to_univ_tnx(params_glwe, pol_univ_computed, pol_computed);

	for (uint64_t p = 0; p < NBASE; p++)
	{
		//TODO: adjust to allow L*KAPPA < 64
		cr_assert(eq(u64, pol_univ[p], pol_univ_computed[p]));
	}

	delete_univ_tnx(pol_univ);
	free(pol_computed);
	delete_univ_tnx(pol_univ_computed);
	delete_glwe_params(params_glwe);
}

//! BIV POLY PART (begin)

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_biv_coef_number returns NBASE * LBASE
	cr_assert(eq(i64, poly_biv_coef_number(params_glwe), NBASE * LBASE,
	             "poly_biv_coef_number failed: got %" PRId64 ", expected %" PRId64, poly_biv_coef_number(params_glwe),
	             NBASE * LBASE));

	// Clean up
	delete_glwe_params(params_glwe);
}

Test(normal_random_biv_poly, does_not_crash)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyBiv* pol = new_biv_poly(params_glwe);

	// Draw normaly pol in Zn[X,Y]
	int status = normal_random_biv_poly(params_glwe, pol);

	// Asserts normal_random_biv_poly passed
	cr_assert(eq(int, status, 0, "normal_random_biv failed."));

	free(pol);
	delete_glwe_params(params_glwe);
}

Test(normal_random_biv_poly, output_is_normalized)
{
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyBiv* pol = new_biv_poly(params_glwe);

	// Draw normaly pol in Zn[X,Y]
	normal_random_biv_poly(params_glwe, pol);

	// Asserts pol is normalized.
	// I.e. that each coefficient is between -2^(KAPPABASE-1) and 2^(KAPPABASE-1)
	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < params_glwe->nn; p++)
		{
			cr_assert(le(i64, pol[(i - 1) * NBASE + p], (1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
			cr_assert(ge(i64, pol[(i - 1) * NBASE + p], -(1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
		}
	}

	free(pol);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
Test(add_biv_poly, basic)
{
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyBiv* pol_lhs      = new_biv_poly(params_glwe);
	PolyBiv* pol_rhs      = new_biv_poly(params_glwe);
	PolyBiv* sum_computed = new_biv_poly(params_glwe);

	normal_random_biv_poly(params_glwe, pol_lhs);
	normal_random_biv_poly(params_glwe, pol_rhs);

	add_biv_poly(module, params_glwe, sum_computed, pol_lhs, pol_rhs);

	// Asserts sum_computed = pol_lhs + pol_rhs
	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < NBASE; p++)
		{
			int64_t idx = p + (i - 1) * params_glwe->nn;
			cr_assert(eq(dbl, sum_computed[idx], pol_lhs[idx] + pol_rhs[idx]),
			          "add_biv_poly mismatch at index %" PRId64 ": %" PRId64 " + %" PRId64 " = %" PRId64
			          ", got %" PRId64,
			          (long long)idx, pol_lhs[idx], pol_rhs[idx], pol_lhs[idx] + pol_rhs[idx], sum_computed[idx]);
		}
	}

	free(pol_lhs);
	free(pol_rhs);
	free(sum_computed);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

//! BIV POLY IN DFT PART (begin)

/**
 * @brief Test normal_random_biv_poly_dft
 *
 *  TODO: what is it even doing? It generates something, normalizes and then checks that it is normalized?
 */
Test(normal_random_biv_poly_dft, unknown_test)
{
	// Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Variables
	PolyBivDFT* pol_dft     = new_biv_poly_dft(params_glwe);
	PolyBiv* pol            = new_biv_poly(params_glwe);
	PolyBiv* pol_normalized = new_biv_poly(params_glwe);

	normal_random_biv_poly_dft(module, params_glwe, pol_dft);

	pvda_vec_znx_idft(module, pol, LBASE, pol_dft, LBASE);

	pvda_vec_znx_normalize_base2k(module, KAPPABASE, pol_normalized, LBASE, NBASE, pol, LBASE, NBASE);

	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < params_glwe->nn; p++)
		{
			cr_assert(le(i64, pol_normalized[(i - 1) * NBASE + p], (1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
			cr_assert(ge(i64, pol_normalized[(i - 1) * NBASE + p], -(1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
		}
	}

	free(pol_dft);
	free(pol);
	free(pol_normalized);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}
