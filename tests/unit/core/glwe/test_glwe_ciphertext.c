#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "core/glwe/glwe_ciphertext.h"
#include "glwe_arithmetic.h"
#include "glwe_params.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

#define NLIMBSBASE (params_glwe->k * glwe_params_l_a(params_glwe) + glwe_params_l_b(params_glwe))
//! COMMON PART (begin)

/**
 * @brief Tests whether glwe_size computes the right size of a bivGLWE ciphertext.
 */
PvdaParamTest(glwe_params_n_limbs, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	cr_assert(eq(i64, glwe_params_n_limbs(params_glwe), NLIMBSBASE));

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether glwe_bytes computes the right number of bytes in a bivGLWE ciphertext.
 */
PvdaParamTest(glwe_bytes, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	cr_assert(eq(i64, glwe_params_bytes(params_glwe), NLIMBSBASE * params_glwe->nn * sizeof(int64_t)));

	DELETE_PVDA_PARAMS_GLWE;
}

//! bivGLWE PART (begin)

/**
 * @brief Tests whether glwe_bytes computes the right number of coefficient in a bivGLWE ciphertext.
 */
PvdaParamTest(glwe_coef_number, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Asserts glwe_coef_number returns NLIMBSBASE * params_glwe->nn
	cr_assert(eq(i64, glwe_coef_number(params_glwe), NLIMBSBASE * params_glwe->nn));

	// Clean up
	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether new_glwe returns a non-NULL pointer.
 */
PvdaParamTest(new_glwe, basic, default_params_fn)
{
	// Parameters
	INIT_PVDA_PARAMS_GLWE(param);
	GLWECiphertext* glwe = new_glwe(params_glwe);

	// Asserts new_glwe allocates a non-NULL glwe
	cr_assert(eq(int, (glwe != NULL) && (glwe->vec != NULL), 1));

	// Clean up
	delete_glwe(glwe);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether add_glwe adds two bivGLWE ciphertexts.
 *
 */
PvdaParamTest(add_glwe, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Variables
	GLWECiphertext* glwe_lhs     = new_glwe(params_glwe);
	GLWECiphertext* glwe_rhs     = new_glwe(params_glwe);
	GLWECiphertext* sum_computed = new_glwe(params_glwe);

	// Draws in Zn[X] the bivGLWE's bivariate elements
	uniform_random_vec(params_glwe->nn, glwe_lhs->vec, glwe_params_n_limbs(params_glwe), params_glwe->nn,
	                   params_glwe->kappa - 1);
	uniform_random_vec(params_glwe->nn, glwe_rhs->vec, glwe_params_n_limbs(params_glwe), params_glwe->nn,
	                   params_glwe->kappa - 1);

	// Computes glwe_lhs + glwe_rhs
	add_glwe(module, sum_computed, glwe_lhs, glwe_rhs);

	// Asserts sum_computed = glwe_lhs + glwe_rhs
	for (uint64_t t = 0; t < glwe_coef_number(params_glwe); t++)
		cr_assert(eq(i64, sum_computed->vec[t], glwe_lhs->vec[t] + glwe_rhs->vec[t]));

	// Clean up
	delete_glwe(glwe_lhs);
	delete_glwe(glwe_rhs);
	delete_glwe(sum_computed);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether sub_glwe subtracts two bivGLWE ciphertexts.
 *
 */
PvdaParamTest(sub_glwe, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Variables
	GLWECiphertext* glwe_lhs     = new_glwe(params_glwe);
	GLWECiphertext* glwe_rhs     = new_glwe(params_glwe);
	GLWECiphertext* sub_observed = new_glwe(params_glwe);

	// Draws in Zn[X] the bivGLWE's bivariate elements
	uniform_random_vec(params_glwe->nn, glwe_lhs->vec, glwe_params_n_limbs(params_glwe), params_glwe->nn,
	                   params_glwe->kappa - 1);
	uniform_random_vec(params_glwe->nn, glwe_rhs->vec, glwe_params_n_limbs(params_glwe), params_glwe->nn,
	                   params_glwe->kappa - 1);

	// Computes glwe_lhs - glwe_rhs
	sub_glwe(module, sub_observed, glwe_lhs, glwe_rhs);

	// Asserts sum_computed = glwe_lhs - glwe_rhs
	for (uint64_t t = 0; t < glwe_coef_number(params_glwe); t++)
		cr_assert(eq(i64, sub_observed->vec[t], glwe_lhs->vec[t] - glwe_rhs->vec[t]));

	// Clean up
	delete_glwe(glwe_lhs);
	delete_glwe(glwe_rhs);
	delete_glwe(sub_observed);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether const_mult_glwe multiply a bivGLWE ciphertext by a ZnX polynomial.
 */
PvdaParamTest(const_mult_glwe, without_normalization, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	// Variables
	GLWECiphertext* prod_computed = new_glwe(params_glwe);
	GLWECiphertext* glwe          = new_glwe(params_glwe);
	PolyUniv* u                   = new_univ(params_glwe);
	PolyUnivDFT* u_dft            = new_univ_dft(module);
	PolyUniv* prod_expected       = new_univ(params_glwe);

	// Draws uniformly the bivGLWE ciphertext and the ZnX polynomial
	uniform_random_vec(params_glwe->nn, glwe->vec, glwe_params_n_limbs(params_glwe), params_glwe->nn,
	                   params_glwe->kappa - 1);

	// Draws in Zn[X] the polynomial u
	uniform_random_pol_znx(u, params_glwe->nn, params_glwe->kappa - 1);

	// Computes u in the DFT domain
	univ_coefs_to_dft(module, u_dft, u);

	// Computes u * glwe
	const_mult_glwe(module, prod_computed, u_dft, glwe);

	// Asserts prod_computed = u * glwe
	for (uint64_t ij = 0; ij < glwe_params_n_limbs(params_glwe); ++ij)
	{
		uint64_t j        = (ij % (params_glwe->k + 1));
		uint64_t i        = ij / (params_glwe->k + 1) + 1;
		PolyUniv* glwe_ij = glwe->vec + (i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn;
		pvda_znx_small_product(module, prod_expected, u, glwe_ij);
		for (uint64_t p = 0; p < params_glwe->nn; p++)
		{
			cr_assert(eq(i64,
			             prod_computed->vec[(i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn + p],
			             prod_expected[p]));
		}
	}

	// Clean up
	delete_univ(u);
	delete_univ_dft(u_dft);
	delete_glwe(glwe);
	delete_glwe(prod_computed);
	delete_univ(prod_expected);

	DELETE_PVDA_PARAMS_GLWE;
}

//! bivGLWE IN DFT PART (begin)

/**
 * @brief Tests whether glwe_coef_number_dft computes the right number of coefficient in a bivGLWE ciphertext in the DFT domain.
 */
PvdaParamTest(glwe_coef_number_dft, basic, default_params_fn)
{
	// Parameters
	INIT_PVDA_PARAMS_GLWE(param);

	// Asserts glwe_coef_number_dft returns NLIMBSBASE * params_glwe->nn / 2
	cr_assert(eq(i64, glwe_coef_number_dft(params_glwe), NLIMBSBASE * params_glwe->nn / 2));

	// Clean up
	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether new_glwe_dft returns a non-NULL pointer.
 */
PvdaParamTest(new_glwe_dft, basic, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	GLWECiphertextDFT* glwe = new_glwe_dft(params_glwe);
	cr_assert(eq(int, (glwe != NULL) && (glwe->vec != NULL), 1));

	delete_glwe_dft(glwe);

	DELETE_PVDA_PARAMS_GLWE;
}

PvdaParamTest(add_glwe_dft, basic, default_params_fn)
{
	// Parameters
	INIT_PVDA_PARAMS_GLWE(param);

	// Variables
	GLWECiphertextDFT* glwe_lhs_dft     = new_glwe_dft(params_glwe);
	GLWECiphertextDFT* glwe_rhs_dft     = new_glwe_dft(params_glwe);
	GLWECiphertextDFT* sum_computed_dft = new_glwe_dft(params_glwe);

	// Draws uniformly the bivGLWE ciphertexts glwe_lhs_dft and glwe_rhs_dft in the DFT domain
	uniform_random_vec_znx_dft(module, glwe_lhs_dft->vec, glwe_params_n_limbs(params_glwe), params_glwe->kappa - 1);
	uniform_random_vec_znx_dft(module, glwe_rhs_dft->vec, glwe_params_n_limbs(params_glwe), params_glwe->kappa - 1);

	// Computes glwe_lhs_dft + glwe_rhs_dft
	add_glwe_dft(sum_computed_dft, glwe_lhs_dft, glwe_rhs_dft);

	// Asserts sum_computed_dft = glwe_lhs_dft + glwe_rhs_dft
	for (uint64_t p = 0; p < glwe_params_n_limbs(params_glwe) * params_glwe->nn; p++)
		cr_assert(eq(sum_computed_dft->vec[p], glwe_lhs_dft->vec[p] + glwe_rhs_dft->vec[p]));

	delete_glwe_dft(glwe_lhs_dft);
	delete_glwe_dft(glwe_rhs_dft);
	delete_glwe_dft(sum_computed_dft);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Tests whether const_mult_glwe_dft multiply a bivGLWE ciphertext by a ZnX polynomial.
 */
PvdaParamTest(const_mult_glwe_dft, without_normalization, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

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
	uniform_random_vec_znx_dft(module, glwe_dft->vec, glwe_params_n_limbs(params_glwe), params_glwe->kappa - 1);

	// Draws uniformly
	uniform_random_pol_znx(u, params_glwe->nn, params_glwe->kappa - 1);

	//! Computation with functions
	// Computes glwe_dft's vec out of the DFT domain
	pvda_vec_znx_idft(module, glwe_ct->vec, glwe_params_n_limbs(params_glwe), glwe_dft->vec,
	                  glwe_params_n_limbs(params_glwe));

	// Computes u in the DFT domain
	univ_coefs_to_dft(module, u_dft, u);

	// Computes DFT(u * glwe)
	const_mult_glwe_dft(module, prod_computed_dft, u_dft, glwe_dft);

	// Computes prod_computed_dft's vec out of the DFT domain
	pvda_vec_znx_idft(module, prod->vec, glwe_params_n_limbs(params_glwe), prod_computed_dft->vec,
	                  glwe_params_n_limbs(params_glwe));

	// Asserts prod_computed_dft = DFT(u * glwe), ie that prod_computed_vec = u * glwe
	for (uint64_t ij = 0; ij < glwe_params_n_limbs(params_glwe); ++ij)
	{
		uint64_t j        = (ij % (params_glwe->k + 1));
		uint64_t i        = ij / (params_glwe->k + 1) + 1;
		PolyUniv* glwe_ij = glwe_ct->vec + (i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn;
		pvda_znx_small_product(module, prod_expected, u, glwe_ij);
		for (uint64_t p = 0; p < params_glwe->nn; p++)
		{
			cr_assert(eq(i64, prod->vec[(i - 1) * (params_glwe->k + 1) * params_glwe->nn + j * params_glwe->nn + p],
			             prod_expected[p]));
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

	DELETE_PVDA_PARAMS_GLWE;
}
