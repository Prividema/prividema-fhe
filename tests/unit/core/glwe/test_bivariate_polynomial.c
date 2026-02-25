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
#define NLIMBSBASE (KBASE+1)*4
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -7

//! COMMON PART (begin)

/**
 * @brief Test poly_biv_size
 */
Test(poly_biv_size, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	cr_assert(eq(i64, poly_biv_size(params), LBASE, "poly_biv_size failed: got %" PRId64 ", expected %" PRId64,
	             poly_biv_size(params), LBASE));

	delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_univ_bytes, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	cr_assert(eq(i64, poly_univ_bytes(params), NBASE * sizeof(int64_t),
	             "poly_univ_bytes failed: got %" PRId64 ", expected %" PRId64, poly_univ_bytes(params),
	             NBASE * sizeof(int64_t)));

	delete_glwe_ct_params(params);
}

/**
 * @brief Test poly_biv_bytes
 */
Test(poly_biv_bytes, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	cr_assert(eq(i64, poly_biv_bytes(params), NBASE * LBASE * sizeof(int64_t),
	             "poly_biv_bytes failed: got %" PRId64 ", expected %" PRId64, poly_biv_bytes(params),
	             NBASE * LBASE * sizeof(int64_t)));

	delete_glwe_ct_params(params);
}

/**
 * @brief Test biv_to_univ
 */
Test(biv_to_univ, test_with_random_biv_generation)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module       = new_module_info_p(NBASE);

	double* pol_univ     = malloc(poly_univ_bytes(params));
	PolyBiv* pol_biv     = new_normal_random_biv_poly(module, params);

	biv_to_univ(params, pol_univ, pol_biv);

	free(pol_univ);
	free(pol_biv);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, one_test)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module       = new_module_info_p(NBASE);

	double* pol_univ     = calloc(poly_univ_bytes(params), 1);
	pol_univ[0]          = 0.0625;

	PolyBiv* pol_biv     = malloc(poly_biv_bytes(params));
	univ_to_biv(params, pol_biv, pol_univ);

	for (int64_t p = 0; p < NBASE; p++) {
		cr_log_info("A %e X^%ld", pol_univ[p], p);
		for (int64_t i = 1; i <= LBASE; i++) cr_log_info("A(XY) %ld Y^%ld", pol_biv[(i - 1) * NBASE + p], i);
	}

	for (int64_t p = 0; p < NBASE; p++) {
		double acc = 0;
		for (int64_t i = 1; i <= LBASE; i++) {
			acc += ldexp((double)pol_biv[(i - 1) * NBASE + p], -i * KAPPABASE);
		}
		cr_log_info("acc %lf pol %lf p %ld", acc, pol_univ[p], p);
		cr_assert(epsilon_eq(dbl, acc - floor(acc) - pol_univ[p] + floor(pol_univ[p]), 0,
		                     ldexp(1.0, -(LBASE - 1) * KAPPABASE)),
		          "acc %lf pol %lf p %ld", acc - floor(acc), -pol_univ[p] + floor(pol_univ[p]), p);
	}
	free(pol_univ);
	free(pol_biv);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module       = new_module_info_p(NBASE);

	double* pol_univ     = malloc(poly_univ_bytes(params));
	new_normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, 1e-2);

	int64_t mask = (1LL << KAPPABASE) - 1;
	for (int64_t p = 0; p < NBASE; p++) {
		cr_log_info("A %e X^%ld", pol_univ[p], p);
		for (int64_t i = 1; i <= LBASE; i++) {
			// cr_log_info("A(XY) %e Y^%ld", ldexp(pol_univ[p], i*kapPABASE), i) ;
			// cr_log_info("A(XY) %ld Y^%ld", (int64_t) ldexp(pol_univ[p], i*kapPABASE) & mask, i) ;
		}
	}

	PolyBiv* pol_biv = malloc(poly_biv_bytes(params));
	univ_to_biv(params, pol_biv, pol_univ);

	for (int64_t p = 0; p < NBASE; p++) {
		cr_log_info("A %e X^%ld", pol_univ[p], p);
		for (int64_t i = 1; i <= LBASE; i++) cr_log_info("A(XY) %ld Y^%ld", pol_biv[(i - 1) * NBASE + p], i);
	}

	for (int64_t p = 0; p < NBASE; p++) {
		double acc = 0;
		for (int64_t i = 1; i <= LBASE; i++) {
			acc += ldexp((double)pol_biv[(i - 1) * NBASE + p], -i * KAPPABASE);
		}
		cr_log_info("acc %lf pol %lf p %ld", acc, pol_univ[p], p);
		cr_assert(epsilon_eq(dbl, acc - floor(acc) - pol_univ[p] + floor(pol_univ[p]), 0,
		                     ldexp(1.0, -(LBASE - 1) * KAPPABASE)),
		          "acc %lf pol %lf p %ld", acc - floor(acc), -pol_univ[p] + floor(pol_univ[p]), p);
	}
	free(pol_univ);
	free(pol_biv);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, maths_test)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module       = new_module_info_p(NBASE);

	// A univariate polynomial drawn normaly in Rn[X]
	double* pol_input_univ = malloc(poly_univ_bytes(params));
	new_normal_random_vec(NBASE, pol_input_univ, 1, NBASE, 0.0, 1e-2);

	PolyBiv* pol_input = malloc(poly_biv_bytes(params));
	univ_to_biv(params, pol_input, pol_input_univ);

	double* pol_computed_univ = calloc(poly_univ_bytes(params), 1);
	biv_to_univ(params, pol_computed_univ, pol_input);

	for (int64_t p = 0; p < NBASE; p++) {
		cr_assert(
		    epsilon_eq(
		        dbl, pol_input_univ[p] - floor(pol_input_univ[p]) - pol_computed_univ[p] + floor(pol_computed_univ[p]),
		        0, ldexp(1.0, -(LBASE - 1) * KAPPABASE)),
		    "pol_input_univ[%ld] %lf pol_computed_univ[%ld] %lf", p, pol_input_univ[p], p, pol_computed_univ[p]);
	}

	free(pol_input_univ);
	free(pol_input);
	free(pol_computed_univ);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

//! BIV POLY PART (begin)

/**
 * @brief Tests if it returns the right size for N in [1,100]
 * @note n_limbs = (k + 1) * l
 */
Test(poly_biv_coef_number, classic_params)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	cr_assert(eq(i64, poly_biv_coef_number(params), NBASE * LBASE,
	             "poly_biv_coef_number failed: got %" PRId64 ", expected %" PRId64, poly_biv_coef_number(params),
	             NBASE * LBASE));

	delete_glwe_ct_params(params);
}

/**
 * @brief Test normal_bivariate_poly
 *
 */
Test(new_normal_random_biv_poly, basic)
{
	MODULE* module       = new_module_info_p(NBASE);
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	PolyBiv* a           = new_normal_random_biv_poly(module, params);

	cr_assert(eq(int, a != NULL, 1, "new_normal_random_biv returned a NULL pointer."));

	free(a);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

/**
 * @brief Test new_normal_random_biv_poly
 *
 */
Test(new_normal_random_biv_poly, is_it_working)
{
	MODULE* module       = new_module_info_p(NBASE);
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	PolyBiv* a           = new_normal_random_biv_poly(module, params);

	for (int64_t i = 1; i <= LBASE; i++) {
		for (int64_t p = 0; p < params->N; p++) {
			cr_assert(le(i64, a[(i - 1) * NBASE + p], (1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
			cr_assert(ge(i64, a[(i - 1) * NBASE + p], -(1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
		}
	}

	free(a);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

/**
 * @brief Test add_biv_poly correctness with random normal polynomials
 */
Test(add_biv_poly, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module       = new_module_info_p(NBASE);

	PolyBiv* a           = new_normal_random_biv_poly(module, params);
	PolyBiv* b           = new_normal_random_biv_poly(module, params);
	PolyBiv* res         = malloc(poly_biv_bytes(params));

	add_biv_poly(params, res, params->N, a, params->N, b, params->N);

	for (int64_t i = 1; i <= LBASE; i++) {
		for (int64_t p = 0; p < NBASE; p++) {
			int64_t idx = p + (i - 1) * params->N;
			cr_assert(eq(dbl, res[idx], a[idx] + b[idx]),
			          "add_biv_poly mismatch at index %" PRId64 ": %" PRId64 " + %" PRId64 " = %" PRId64
			          ", got %" PRId64,
			          (long long)idx, a[idx], b[idx], a[idx] + b[idx], res[idx]);
		}
	}

	free(a);
	free(b);
	free(res);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}

//! BIV POLY IN DFT PART (begin)

/**
 * @brief Test add_biv_poly_dft correctness with random normal DFT polynomials
 */
Test(add_biv_poly_dft, basic)
{
	GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
	MODULE* module       = new_module_info_p(NBASE);

	PolyBivDFT* a        = new_normal_random_biv_poly_dft(module, params);
	PolyBivDFT* b        = new_normal_random_biv_poly_dft(module, params);
	PolyBivDFT* res      = malloc(poly_biv_bytes(params));

	add_biv_poly_dft(params, res, params->N, a, params->N, b, params->N);

	for (int64_t i = 1; i <= LBASE; i++) {
		for (int64_t p = 0; p < NBASE; p++) {
			int64_t idx = p + (i - 1) * params->N;
			cr_assert(epsilon_eq(dbl, res[idx], a[idx] + b[idx], 1e-9),
			          "add_biv_poly_dft mismatch at index %" PRId64 ": %f + %f = %f, got %f", (long long)idx, a[idx],
			          b[idx], a[idx] + b[idx], res[idx]);
		}
	}

	free(a);
	free(b);
	free(res);
	delete_module_info_p(module);
	delete_glwe_ct_params(params);
}

/**
 * @brief Test new_normal_random_biv_poly_dft
 *
 */
Test(new_normal_random_biv_poly_dft, is_it_working)
{
	MODULE* module        = new_module_info_p(NBASE);
	GLWECtParams* params  = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));

	PolyBivDFT* a_dft     = new_normal_random_biv_poly_dft(module, params);
	PolyBiv* a            = malloc(poly_biv_bytes(params));
	PolyBiv* a_normalized = malloc(poly_biv_bytes(params));

	vec_znx_idft_p(module, a, LBASE, a_dft, LBASE);
	vec_znx_normalize_base2k_p(module, KAPPABASE, a_normalized, LBASE, NBASE, a, LBASE, NBASE);

	for (int64_t i = 1; i <= LBASE; i++) {
		for (int64_t p = 0; p < params->N; p++) {
			cr_assert(le(i64, a_normalized[(i - 1) * NBASE + p], (1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is greater than 2^(kappa-1)."));
			cr_assert(ge(i64, a_normalized[(i - 1) * NBASE + p], -(1 << (KAPPABASE - 1)),
			             "The coefficient of a(X^p, Y^i) is smaller than -2^(kappa-1)."));
		}
	}
	cr_assert(1);

	free(a_dft);
	free(a);
	free(a_normalized);
	delete_glwe_ct_params(params);
	delete_module_info_p(module);
}