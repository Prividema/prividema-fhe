#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <inttypes.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "bivariate_polynomial.h"
#include "common/rng.h"
#include "common/spqlios_alias.h"
#include "core/glwe/glwe.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 4
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -7

//! COMMON PART (begin)

/**
 * @brief Test poly_biv_size
 */
Test(poly_biv_size, basic)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_biv_size returns LBASE
	cr_assert(eq(i64, poly_biv_size(params_glwe), LBASE, "poly_biv_size failed: got %" PRId64 ", expected %" PRId64,
	             poly_biv_size(params_glwe), LBASE));

	// Clean up
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_univ_bytes, basic)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_univ_bytes returns NBASE * sizeof(int64_t)
	cr_assert(eq(i64, poly_univ_bytes(params_glwe), NBASE * sizeof(int64_t),
	             "poly_univ_bytes failed: got %" PRId64 ", expected %" PRId64, poly_univ_bytes(params_glwe),
	             NBASE * sizeof(int64_t)));

	// Clean up
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_biv_bytes, basic)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_biv_bytes returns NBASE * LBASE * sizeof(int64_t)
	cr_assert(eq(i64, poly_biv_bytes(params_glwe), NBASE * LBASE * sizeof(int64_t),
	             "poly_biv_bytes failed: got %" PRId64 ", expected %" PRId64, poly_biv_bytes(params_glwe),
	             NBASE * LBASE * sizeof(int64_t)));

	// Clean up
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Test biv_to_univ
 */
Test(biv_to_univ, runs)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module            = pvda_new_module_info(NBASE);

	// Variables
	PolyBiv* pol                   = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* pol_univ_computed = malloc(poly_univ_bytes(params_glwe));

	// Define pol_biv = Y
	pol[0] = 1;

	// Computes the bivariate polynomial in the Tn[X]
	biv_to_univ(params_glwe, pol_univ_computed, pol);

	// Asserts pol_univ_computed(X) = 2^{-KAPPABASE}
	cr_assert(eq(dbl, pol_univ_computed[0], ldexp(1.0, -KAPPABASE)));

	for (uint64_t p = 1; p < NBASE; p++)
		cr_assert(eq(dbl, pol_univ_computed[p], 0), "pol_univ_computed[%ld] = %ld ", p, pol_univ_computed[p]);

	// Clean up
	free(pol);
	free(pol_univ_computed);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, one_test)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module            = pvda_new_module_info(NBASE);

	// Variables
	double* pol_univ      = calloc(poly_univ_bytes(params_glwe), 1);
	PolyBiv* pol_computed = malloc(poly_biv_bytes(params_glwe));

	// Define pol_univ(X) = 2^{-KAPPABASE} , i.e. in Zn[X,Y] pol(X,Y) = Y
	pol_univ[0] = ldexp(1.0, -KAPPABASE);

	// Compute pol_univ's base-2KAPPABASE normalized decomposition
	univ_to_biv(params_glwe, pol_computed, pol_univ);

	// Asserts pol_computed = Y
	cr_assert(eq(int, pol_computed[0], 1), "pol_computed[%ld, %ld] = %ld ", 0, 1, pol_computed[0]);

	for (uint64_t p = 1; p < NBASE; p++)
		cr_assert(eq(int, pol_computed[p], 0), "pol_computed[%ld, %ld] = %ld ", p, 1, pol_computed[p]);

	for (uint64_t i = 2; i <= LBASE; i++)
		for (uint64_t p = 0; p < NBASE; p++)
			cr_assert(eq(int, pol_computed[(i - 1) * NBASE + p], 0), "pol_computed[%ld, %ld] = %ld ", p, i,
			          pol_computed[(i - 1) * NBASE + p]);

	// Clean up
	free(pol_univ);
	free(pol_computed);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, basic)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module            = new_module_info(NBASE, FFT64);

	// Variables
	double* pol_univ      = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* pol_computed = malloc(poly_biv_bytes(params_glwe));

	// Draws normaly pol_univ in Rn[X]
	normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, 1e-2);

	// Computes pol_univ's base-2KAPPABASE normalized decomposition
	univ_to_biv(params_glwe, pol_computed, pol_univ);

	// pol_computed_p = pol_computed(X, Y = 2^{-KAPPABASE})[p]
	double pol_computed_p = 0;

	// Asserts pol_computed, in Rn[X] (with Y = 2^{-KAPPABASE}), is equal to pol_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{
		pol_computed_p = 0;
		for (uint64_t i = 1; i <= LBASE; i++)
			pol_computed_p += ldexp((double)pol_computed[(i - 1) * NBASE + p], -i * KAPPABASE);

		cr_assert(epsilon_eq(dbl, pol_computed_p - floor(pol_computed_p) - pol_univ[p] + floor(pol_univ[p]), 0,
		                     ldexp(1.0, -(LBASE - 1) * KAPPABASE)),
		          "pol_univ(X = %ld) = %lf and pol(X = %ld) %lf", p, pol_computed_p - floor(pol_computed_p),
		          -pol_univ[p] + floor(pol_univ[p]), p);
	}

	// Clean up
	free(pol_univ);
	free(pol_computed);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, maths_test)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module            = pvda_new_module_info(NBASE);

	// Variables
	double* pol_univ          = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* pol_computed     = malloc(poly_biv_bytes(params_glwe));
	double* pol_univ_computed = calloc(poly_univ_bytes(params_glwe), 1);

	// Draw pol_univ normaly in Rn[X]
	normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, 1e-2);

	// Computes pol_univ's base-2KAPPABASE decomposition
	univ_to_biv(params_glwe, pol_computed, pol_univ);

	// Computes pol in Rn[X] with Y = 2^{-KAPPABASE}
	biv_to_univ(params_glwe, pol_univ_computed, pol_computed);

	// Asserts pol_univ_computed(X) = pol_univ(X)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		cr_assert(epsilon_eq(dbl, pol_univ[p] - floor(pol_univ[p]) - pol_univ_computed[p] + floor(pol_univ_computed[p]),
		                     0, ldexp(1.0, -(LBASE - 1) * KAPPABASE)),
		          "pol_univ[%ld] = %lf and pol_univ_computed[%ld] = %lf", p, pol_univ[p], p, pol_univ_computed[p]);
	}

	// Clean up
	free(pol_univ);
	free(pol_computed);
	free(pol_univ_computed);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

//! BIV POLY PART (begin)

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params)
{
	// Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Asserts poly_biv_coef_number returns NBASE * LBASE
	cr_assert(eq(i64, poly_biv_coef_number(params_glwe), NBASE * LBASE,
	             "poly_biv_coef_number failed: got %" PRId64 ", expected %" PRId64, poly_biv_coef_number(params_glwe),
	             NBASE * LBASE));

	// Clean up
	delete_glwe_ct_params(params_glwe);
}

/**
 * @brief Test normal_bivariate_poly
 *
 */
Test(normal_random_biv_poly, basic)
{
	// Parameters
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Variables
	PolyBiv* pol = malloc(poly_biv_bytes(params_glwe));
	;

	// Draw normaly pol in Zn[X,Y]
	int status = normal_random_biv_poly(params_glwe, pol);

	// Asserts normal_random_biv_poly passed
	cr_assert(eq(int, status, 0, "normal_random_biv failed."));

	// Clean up
	free(pol);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test normal_random_biv_poly
 *
 */
Test(normal_random_biv_poly, is_it_working)
{
	// Parameters
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Variables
	PolyBiv* pol = malloc(poly_biv_bytes(params_glwe));
	;

	// Draw normaly pol in Zn[X,Y]
	normal_random_biv_poly(params_glwe, pol);

	// Asserts pol is normalized.
	// I.e. that each coefficient is between -2^(KAPPABASE-1) and 2^(KAPPABASE-1)
	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < params_glwe->N; p++)
		{
			cr_assert(le(i64, pol[(i - 1) * NBASE + p], (1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
			cr_assert(ge(i64, pol[(i - 1) * NBASE + p], -(1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
		}
	}

	// Clean up
	free(pol);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
Test(add_biv_poly, basic)
{
	// Parameters
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Variables
	PolyBiv* pol_lhs = malloc(poly_biv_bytes(params_glwe));
	;
	PolyBiv* pol_rhs = malloc(poly_biv_bytes(params_glwe));
	;
	PolyBiv* sum_computed = malloc(poly_biv_bytes(params_glwe));

	// Draw normaly pol_lhs and pol_rhs in Zn[X,Y]
	normal_random_biv_poly(params_glwe, pol_lhs);
	normal_random_biv_poly(params_glwe, pol_rhs);

	// Computes pol_lhs + pol_rhs
	add_biv_poly(params_glwe, sum_computed, params_glwe->N, pol_lhs, params_glwe->N, pol_rhs, params_glwe->N);

	// Asserts sum_computed = pol_lhs + pol_rhs
	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < NBASE; p++)
		{
			int64_t idx = p + (i - 1) * params_glwe->N;
			cr_assert(eq(dbl, sum_computed[idx], pol_lhs[idx] + pol_rhs[idx]),
			          "add_biv_poly mismatch at index %" PRId64 ": %" PRId64 " + %" PRId64 " = %" PRId64
			          ", got %" PRId64,
			          (long long)idx, pol_lhs[idx], pol_rhs[idx], pol_lhs[idx] + pol_rhs[idx], sum_computed[idx]);
		}
	}

	// Clean up
	free(pol_lhs);
	free(pol_rhs);
	free(sum_computed);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

//! BIV POLY IN DFT PART (begin)

/**
 * @brief Test normal_random_biv_poly_dft
 *
 */
Test(normal_random_biv_poly_dft, is_it_working)
{
	// Parameters
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Variables
	PolyBivDFT* pol_dft = malloc(poly_biv_bytes(params_glwe));
	;
	PolyBiv* pol            = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* pol_normalized = malloc(poly_biv_bytes(params_glwe));

	// Draw pol_dft normaly in the DFT domain
	normal_random_biv_poly_dft(module, params_glwe, pol_dft);

	// Computes pol_dft out of the DFT domain
	pvda_vec_znx_idft(module, pol, LBASE, pol_dft, LBASE);

	// Normalize pol
	pvda_vec_znx_normalize_base2k(module, KAPPABASE, pol_normalized, LBASE, NBASE, pol, LBASE, NBASE);

	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < params_glwe->N; p++)
		{
			cr_assert(le(i64, pol_normalized[(i - 1) * NBASE + p], (1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
			cr_assert(ge(i64, pol_normalized[(i - 1) * NBASE + p], -(1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
		}
	}
	cr_assert(1);

	free(pol_dft);
	free(pol);
	free(pol_normalized);
	delete_glwe_ct_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test add_biv_poly_dft correctness with random normal DFT polynomials
 */
Test(add_biv_poly_dft, basic)
{
	// Parameters
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	// Variables
	PolyBivDFT* pol_lhs_dft = malloc(poly_biv_bytes(params_glwe));
	;
	PolyBivDFT* pol_rhs_dft = malloc(poly_biv_bytes(params_glwe));
	;
	PolyBivDFT* sum_computed_dft = malloc(poly_biv_bytes(params_glwe));

	// Draw normaly pol_lhs and pol_rhs in Zn[X,Y]
	normal_random_biv_poly_dft(module, params_glwe, pol_lhs_dft);
	normal_random_biv_poly_dft(module, params_glwe, pol_rhs_dft);

	add_biv_poly_dft(params_glwe, sum_computed_dft, params_glwe->N, pol_lhs_dft, params_glwe->N, pol_rhs_dft,
	                 params_glwe->N);

	for (uint64_t i = 1; i <= LBASE; i++)
	{
		for (uint64_t p = 0; p < NBASE; p++)
		{
			int64_t idx = p + (i - 1) * params_glwe->N;
			cr_assert(epsilon_eq(dbl, sum_computed_dft[idx], pol_lhs_dft[idx] + pol_rhs_dft[idx], 1e-9),
			          "add_biv_poly_dft mismatch at index %" PRId64 ": %f + %f = %f, got %f", (long long)idx,
			          pol_lhs_dft[idx], pol_rhs_dft[idx], pol_lhs_dft[idx] + pol_rhs_dft[idx], sum_computed_dft[idx]);
		}
	}

	free(pol_lhs_dft);
	free(pol_rhs_dft);
	free(sum_computed_dft);
	pvda_delete_module_info(module);
	delete_glwe_ct_params(params_glwe);
}
