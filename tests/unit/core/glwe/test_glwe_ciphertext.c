#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "core/glwe/glwe_ciphertext.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "univariate_polynomial.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 2
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  1e-7

//! COMMON PART (begin)

/**
 * @brief Tests whether glwe_size computes the right size of a bivGLWE ciphertext.
 */
Test(glwe_size, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	// Asserts glwe_size returns NLIMBSBASE
	cr_assert(eq(i64, glwe_total_nlimbs(params_glwe), NLIMBSBASE));

	// Clean up
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether glwe_bytes computes the right number of bytes in a bivGLWE ciphertext.
 */
Test(glwe_bytes, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	// Asserts glwe_bytes returns NLIMBSBASE * NBASE * sizeof(int64_t)
	cr_assert(eq(i64, glwe_bytes(params_glwe), NLIMBSBASE * NBASE * sizeof(int64_t)));

	// Clean up
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether mult_vec_znx_dft multiply correctly two Zn[X] polynomials pol_lhs and b. Ie res = pol_lhs * b.
 */
Test(mult_vec_znx_dft, size_equal_one)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module          = pvda_new_module_info(NBASE);

	// Variables
	// TODO: WONTFIX
	int64_t* prod_computed    = calloc(poly_univ_bytes(params_glwe), 1);
	int64_t* pol_lhs          = calloc(poly_univ_bytes(params_glwe), 1);
	int64_t* pol_rhs          = calloc(poly_univ_bytes(params_glwe), 1);
	double* prod_computed_dft = calloc(poly_univ_bytes(params_glwe), 1);
	double* pol_lhs_dft       = calloc(poly_univ_bytes(params_glwe), 1);
	double* pol_rhs_dft       = calloc(poly_univ_bytes(params_glwe), 1);
	int64_t* prod_expected    = calloc(poly_univ_bytes(params_glwe), 1);

	// Draws uniformly in Zn[X]
	uniform_random_pol_znx(pol_lhs, NBASE, 14);
	uniform_random_pol_znx(pol_rhs, NBASE, 14);

	// Computes in the DFT pol_lhs and pol_rhs
	univ_coefs_to_dft(module, pol_lhs_dft, pol_lhs);
	univ_coefs_to_dft(module, pol_rhs_dft, pol_rhs);

	// prod_computed_dft = DFT(pol_lhs ⊙ pol_rhs)
	mult_vec_znx_dft(module, prod_computed_dft, 1, pol_lhs_dft, 1, pol_rhs_dft, 1);

	univ_dft_to_coefs(module, prod_computed, prod_computed_dft);

	// Compare the real coefficient res_p for each p in [0, NBASE -1] with the res_p mult_vec_znx_dft computed
	// coefficient.
	pvda_znx_product(module, prod_expected, pol_lhs, pol_rhs);

	for (uint64_t p = 0; p < NBASE; p++)
	{
		cr_assert(eq(i64, prod_computed[p], prod_expected[p]));
	}

	free(prod_computed);
	free(prod_computed_dft);
	free(prod_expected);
	free(pol_lhs);
	free(pol_lhs_dft);
	free(pol_rhs);
	free(pol_rhs_dft);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether mult_vec_znx_dft correctly multiplies two Zn[X] vectors vec_lhs and vec_rhs component-wise..
 * It draws a random uniform size, ie a random uniform number of Zn[X] polynomials.
 * Ie vec_lhs = (a_i), vec_rhs = (b_i) -> res = vec_lhs ⊙ vec_rhs = (a_i * b_i). Where a_i and b_i are in Zn[X]
 */
Test(mult_vec_znx_dft, random_size)
{
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module          = pvda_new_module_info(NBASE);

	int64_t size = 0;

	while (size <= 0)
	{
		rand_uniform(&size, 8);
	}

	// Variables
	// TODO: WONTFIX
	int64_t* component_wise_mult    = calloc(poly_univ_bytes(params_glwe) * size, 1);
	int64_t* vec_lhs                = calloc(poly_univ_bytes(params_glwe) * size, 1);
	int64_t* vec_rhs                = calloc(poly_univ_bytes(params_glwe) * size, 1);
	double* component_wise_mult_dft = calloc(poly_univ_bytes(params_glwe) * size, 1);
	double* vec_lhs_dft             = calloc(poly_univ_bytes(params_glwe) * size, 1);
	double* vec_rhs_dft             = calloc(poly_univ_bytes(params_glwe) * size, 1);
	int64_t* prod_expected          = calloc(poly_univ_bytes(params_glwe), 1);

	// Draws uniformly in (Zn[X])^size vec_lhs and vec_rhs
	uniform_random_vec(NBASE, vec_lhs, size, NBASE, 14);
	uniform_random_vec(NBASE, vec_rhs, size, NBASE, 14);

	// Computes in the DFT vec_lhs and vec_rhs
	pvda_vec_znx_dft(module, vec_lhs_dft, size, vec_lhs, size, NBASE);
	pvda_vec_znx_dft(module, vec_rhs_dft, size, vec_rhs, size, NBASE);

	// Computes component_wise_mult_dft = DFT(vec_lhs) ⊙ DFT(vec_rhs)
	mult_vec_znx_dft(module, component_wise_mult_dft, size, vec_lhs_dft, size, vec_rhs_dft, size);

	// component_wise_mult = vec_lhs ⊙ vec_rhs
	pvda_vec_znx_idft(module, component_wise_mult, size, component_wise_mult_dft, size);

	// Compare the real coefficient component_wise_mult_p for each p in [0, NBASE -1] with the component_wise_mult_p mult_vec_znx_dft computed
	// coefficient.
	for (uint64_t i = 0; i < size; i++)
	{
		pvda_znx_product(module, prod_expected, vec_lhs + i * NBASE, vec_rhs + i * NBASE);
		for (uint64_t p = 0; p < NBASE; p++)
		{
			cr_assert(eq(i64, component_wise_mult[i * NBASE + p], prod_expected[p]));
		}
	}

	free(component_wise_mult);
	free(component_wise_mult_dft);
	free(vec_lhs);
	free(vec_lhs_dft);
	free(vec_rhs);
	free(vec_rhs_dft);
	free(prod_expected);
	pvda_delete_module_info(module);
	delete_glwe_params(params_glwe);
}

//! bivGLWE PART (begin)

/**
 * @brief Tests whether glwe_bytes computes the right number of coefficient in a bivGLWE ciphertext.
 */
Test(glwe_coef_number, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	// Asserts glwe_coef_number returns NLIMBSBASE * NBASE
	cr_assert(eq(i64, glwe_coef_number(params_glwe), NLIMBSBASE * NBASE));

	// Clean up
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether new_glwe returns a non-NULL pointer.
 */
Test(new_glwe, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GLWECiphertext* glwe    = new_glwe(params_glwe);

	// Asserts new_glwe allocates a non-NULL glwe
	cr_assert(eq(int, (glwe != NULL) && (glwe->vec != NULL), 1));

	// Clean up
	delete_glwe(glwe);
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether add_glwe adds two bivGLWE ciphertexts.
 */
Test(add_glwe, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	// Variables
	GLWECiphertext* glwe_lhs     = new_glwe(params_glwe);
	GLWECiphertext* glwe_rhs     = new_glwe(params_glwe);
	GLWECiphertext* sum_computed = new_glwe(params_glwe);

	// Draws in Zn[X] the bivGLWE's bivariate elements
	uniform_random_vec(NBASE, glwe_lhs->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);
	uniform_random_vec(NBASE, glwe_rhs->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);

	// Computes glwe_lhs + glwe_rhs
	add_glwe(sum_computed, glwe_lhs, glwe_rhs);

	// Asserts sum_computed = glwe_lhs + glwe_rhs
	for (uint64_t t = 0; t < glwe_coef_number(params_glwe); t++)
		cr_assert(eq(sum_computed->vec[t], glwe_lhs->vec[t] + glwe_rhs->vec[t]));

	// Clean up
	delete_glwe(glwe_lhs);
	delete_glwe(glwe_rhs);
	delete_glwe(sum_computed);
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether const_mult_glwe multiply a bivGLWE ciphertext by a ZnX polynomial.
 */
Test(const_mult_glwe, without_normalization)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module          = pvda_new_module_info(NBASE);

	// Variables
	GLWECiphertext* prod_computed = new_glwe(params_glwe);
	GLWECiphertext* glwe          = new_glwe(params_glwe);
	PolyUniv* u                   = new_univ(params_glwe);
	PolyUnivDFT* u_dft            = new_univ_dft(module);
	PolyUniv* prod_expected       = new_univ(params_glwe);

	// Draws uniformly the bivGLWE ciphertext and the ZnX polynomial
	uniform_random_vec(NBASE, glwe->vec, params_glwe->n_limbs, NBASE, KAPPABASE - 1);

	// Draws in Zn[X] the polynomial u
	uniform_random_pol_znx(u, NBASE, KAPPABASE - 1);

	// Computes u in the DFT domain
	univ_coefs_to_dft(module, u_dft, u);

	// Computes u * glwe
	const_mult_glwe(module, prod_computed, u_dft, glwe);

	// Asserts prod_computed = u * glwe
	for (uint64_t i = 1; i <= LBASE; i++)
		for (uint64_t j = 0; j < KBASE + 1; j++)
		{
			PolyUniv* glwe_ij = glwe->vec + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;
			pvda_znx_product(module, prod_expected, u, glwe_ij);
			for (uint64_t p = 0; p < NBASE; p++)
			{
				cr_assert(eq(i64, prod_computed->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p], prod_expected[p]));
			}
		}

	// Clean up
	delete_univ(u);
	delete_univ_dft(u_dft);
	pvda_delete_module_info(module);
	delete_glwe(glwe);
	delete_glwe(prod_computed);
	delete_univ(prod_expected);
	delete_glwe_params(params_glwe);
}

//! bivGLWE IN DFT PART (begin)

/**
 * @brief Tests whether glwe_coef_number_dft computes the right number of coefficient in a bivGLWE ciphertext in the DFT domain.
 */
Test(glwe_coef_number_dft, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);

	// Asserts glwe_coef_number_dft returns NLIMBSBASE * NBASE / 2
	cr_assert(eq(i64, glwe_coef_number_dft(params_glwe), NLIMBSBASE * NBASE / 2));

	// Clean up
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether new_glwe_dft returns a non-NULL pointer.
 */
Test(new_glwe_dft, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	GLWECiphertextDFT* glwe = new_glwe_dft(params_glwe);

	// Asserts new_glwe_dft allocates a non-NULL glwe
	cr_assert(eq(int, (glwe != NULL) && (glwe->vec != NULL), 1));

	// Clean up
	delete_glwe_dft(glwe);
	delete_glwe_params(params_glwe);
}

Test(add_glwe_dft, basic)
{
	// Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module          = pvda_new_module_info(NBASE);

	// Variables
	GLWECiphertextDFT* glwe_lhs_dft     = new_glwe_dft(params_glwe);
	GLWECiphertextDFT* glwe_rhs_dft     = new_glwe_dft(params_glwe);
	GLWECiphertextDFT* sum_computed_dft = new_glwe_dft(params_glwe);

	// Draws uniformly the bivGLWE ciphertexts glwe_lhs_dft and glwe_rhs_dft in the DFT domain
	uniform_random_vec_znx_dft(module, glwe_lhs_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);
	uniform_random_vec_znx_dft(module, glwe_rhs_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);

	// Computes glwe_lhs_dft + glwe_rhs_dft
	add_glwe_dft(sum_computed_dft, glwe_lhs_dft, glwe_rhs_dft);

	// Asserts sum_computed_dft = glwe_lhs_dft + glwe_rhs_dft
	for (uint64_t i = 1; i < LBASE; i++)
		for (uint64_t j = 0; j < KBASE + 1; j++)
			for (uint64_t p = 0; p < NBASE; p++)
				cr_assert(eq(sum_computed_dft->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p],
				             glwe_lhs_dft->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p] +
				                 glwe_rhs_dft->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p]));

	// Clean up
	pvda_delete_module_info(module);
	delete_glwe_dft(glwe_lhs_dft);
	delete_glwe_dft(glwe_rhs_dft);
	delete_glwe_dft(sum_computed_dft);
	delete_glwe_params(params_glwe);
}

/**
 * @brief Tests whether const_mult_glwe_dft multiply a bivGLWE ciphertext by a ZnX polynomial.
 */
Test(const_mult_glwe_dft, without_normalization)
{
	//! Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
	MODULE* module          = pvda_new_module_info(NBASE);

	//! Variables
	GLWECiphertextDFT* prod_computed_dft = new_glwe_dft(params_glwe);

	GLWECiphertext* prod        = new_glwe(params_glwe);
	GLWECiphertextDFT* glwe_dft = new_glwe_dft(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);
	PolyUniv* u                 = new_univ(params_glwe);
	PolyUnivDFT* u_dft          = new_univ_dft(module);
	PolyUniv* prod_expected     = new_univ(params_glwe);

	//! Draws input variables
	// Draws uniformly the bivGLWE ciphertext in the DFT domain
	uniform_random_vec_znx_dft(module, glwe_dft->vec, params_glwe->n_limbs, KAPPABASE - 1);

	// Draws uniformly
	uniform_random_pol_znx(u, NBASE, KAPPABASE - 1);

	//! Computation with functions
	// Computes glwe_dft's vec out of the DFT domain
	pvda_vec_znx_idft(module, glwe_ct->vec, glwe_total_nlimbs(params_glwe), glwe_dft->vec,
	                  glwe_total_nlimbs(params_glwe));

	// Computes u in the DFT domain
	univ_coefs_to_dft(module, u_dft, u);

	// Computes DFT(u * glwe)
	const_mult_glwe_dft(module, prod_computed_dft, u_dft, glwe_dft);

	// Computes prod_computed_dft's vec out of the DFT domain
	pvda_vec_znx_idft(module, prod->vec, glwe_total_nlimbs(params_glwe), prod_computed_dft->vec,
	                  glwe_total_nlimbs(params_glwe));

	// Asserts prod_computed_dft = DFT(u * glwe), ie that prod_computed_vec = u * glwe
	for (uint64_t i = 1; i <= LBASE; i++)
		for (uint64_t j = 0; j < KBASE + 1; j++)
		{
			PolyUniv* glwe_ij = glwe_ct->vec + (i - 1) * (KBASE + 1) * NBASE + j * NBASE;
			pvda_znx_product(module, prod_expected, u, glwe_ij);
			for (uint64_t p = 0; p < NBASE; p++)
			{
				cr_assert(eq(i64, prod->vec[(i - 1) * (KBASE + 1) * NBASE + j * NBASE + p], prod_expected[p]));
			}
		}

	// Clean up
	delete_univ(u);
	delete_univ_dft(u_dft);
	delete_glwe(glwe_ct);
	delete_glwe(prod);
	delete_univ(prod_expected);
	delete_glwe_dft(glwe_dft);
	delete_glwe_dft(prod_computed_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}
