#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "common/spqlios_alias.h"
#include "core/glwe/glwe.h"
#include "rng.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 8
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -(LBASE / 2 + 1) * KAPPABASE

//! GGWS PART (begin)

/**
 * @brief Test glwe_secret_masking. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking, small_error)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	//! Parameters
	MODULE* module            = new_module_info_p(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKeyDFT* sk_dft    = new_glwe_secret_key_dft(NBASE, KBASE);
	GLWECiphertext* glwe_computed        = new_glwe(params_glwe);
	PolyBiv* m                  = malloc(poly_biv_bytes(params_glwe));
	double* m_univ              = calloc(NBASE, sizeof(double));
	PolyBiv* err                = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase              = calloc(NBASE * LBASE, sizeof(int64_t));
	PolyBiv* phase_computed     = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	double* phase_computed_univ = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key_dft(module, sk_dft, 2);

	// The input message uniformly drawn in Zn[X,Y]
	uniform_random_biv_poly(module, params_glwe, m, LBASE / 2);

	// The input error normaly drawn in Zn[X,Y]
	normal_random_biv_poly(module, params_glwe, err);

	//! Computations with functions
	// Computes the message in Tn[X] with the base-2Kappa
	biv_to_univ(params_glwe, m_univ, m);

	// The final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking(module, glwe_computed, sk_dft, phase);

	// The computed phase in Rn[X]
	glwe_secret_demasking(module, phase_computed, sk_dft, glwe_computed);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	//! Asserts phase_computed_univ = m + err
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff_1 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p];
		double diff_2 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3*sigma;

		// The difference should be smaller than : |err_p| Ie 3*sigma, 99.73% of the time
		int cond = (diff_1 <= err_length && diff_1 >= -err_length) || (diff_2 <= err_length && diff_2 >= -err_length);

		if (!cond) big_error_count++;
	}

	// Asserts big_error_count <= 0.27*N
	cr_assert(big_error_count <= (int)(0.27 * NBASE), "The error should be greater than 3*sigma at most %ld times", (int)(0.27*NBASE));

	//! Clean up
	free(m);
	free(m_univ);
	free(err);
	free(phase);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe(glwe_computed);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_ct_params(params_glwe);
	delete_module_info_p(module);
}

/**
 * @brief Test glwe_encrypt_priv. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking, uniform_RnX_message)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	//! Parameters
	MODULE* module            = new_module_info_p(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKeyDFT* sk_dft    = new_glwe_secret_key_dft(NBASE, KBASE);
	GLWECiphertext* glwe_computed   = malloc(glwe_coef_number(params_glwe) * sizeof(int64_t));
	double* m_univ              = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* phase              = calloc(NBASE * LBASE, sizeof(int64_t));
	PolyBiv* m                  = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* err                = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase_computed     = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	double* phase_computed_univ = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key_dft(module, sk_dft, 2);

	// Draws normaly in Rn[X] m_univ
	normal_random_vec(NBASE, m_univ, 1, NBASE, 0.0, 0.1);
	
	// Draws normaly the error
	normal_random_biv_poly(module, params_glwe, err);

	//! Computation with functions
	// Computes m_univ bivariate form
	univ_to_biv(params_glwe, m, m_univ);

	// Computes the final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking(module, glwe_computed, sk_dft, phase);

	// Computes the computed phase in Rn[X]
	glwe_secret_demasking(module, phase_computed, sk_dft, glwe_computed);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	
	//! Asserts phase_computed_univ = m + err
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff_1 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p];
		double diff_2 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma + ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);
		
		// Using the triangle inequality, the difference should be smaller than :
		// |err_p| + |m_p - m_computed_p| <= 3*sigma + 2^(-l*kappa)
		int cond = (diff_1 <= err_length && diff_1 >= -err_length) || (diff_2 <= err_length && diff_2 >= -err_length);

		if (!cond) big_error_count++;
	}

	// Asserts big_error_count <= 0.27*N
	cr_assert(big_error_count <= (int)(0.27 * NBASE), "The error should be greater than 3*sigma at most %ld times", (int)(0.27*NBASE));

	free(m);
	free(m_univ);
	free(err);
	free(phase);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe(glwe_computed);
	delete_module_info_p(module);
	delete_glwe_ct_params(params_glwe);
	delete_glwe_secret_key_dft(sk_dft);
}


//! GLWE IN DFT PART (begin)
/**
 * @brief Test glwe_secret_masking. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking_dft, small_error)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	//! Parameters
	MODULE* module            = new_module_info_p(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKeyDFT* sk_dft      = new_glwe_secret_key_dft(NBASE, KBASE);
	GLWECiphertextDFT* glwe_computed_dft = new_glwe_dft(params_glwe);
	PolyBiv* m                    = malloc(poly_biv_bytes(params_glwe));
	double* m_univ                = calloc(NBASE, sizeof(double));
	PolyBiv* err                  = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase                = calloc(NBASE * LBASE, sizeof(int64_t));
	PolyBivDFT* phase_dft         = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase_computed       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	double* phase_computed_univ   = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key_dft(module, sk_dft, 2);

	// The input message uniformly drawn in Zn[X,Y]
	uniform_random_biv_poly(module, params_glwe, m, LBASE / 2);

	// The input error normaly drawn in Zn[X,Y]
	normal_random_biv_poly(module, params_glwe, err);

	//! Computes with functions
	// Computes the message in Tn[X] with the base-2Kappa
	biv_to_univ(params_glwe, m_univ, m);

	// The final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the phase in the DFT domain
	vec_znx_dft_p(module, phase_dft, LBASE, phase, LBASE, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking_dft(module, glwe_computed_dft, sk_dft, phase_dft);

	// The computed phase in Rn[X]
	glwe_secret_demasking_dft(module, phase_computed, sk_dft, glwe_computed_dft);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Compare both phase in Rn[X]
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff_1 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p];
		double diff_2 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma ;

		// The difference should be smaller than : |err_p| Ie 3*sigma, 99.73% of the time
		int cond = (diff_1 <= err_length && diff_1 >= -err_length) || (diff_2 <= err_length && diff_2 >= -err_length);

		if (!cond) big_error_count++;
	}

	/// Asserts big_error_count <= 0.27*N
	cr_assert(big_error_count <= (int)(0.27 * NBASE), "The error should be greater than 3*sigma at most %ld times", (int)(0.27*NBASE));

	//! Clean up
	free(m);
	free(m_univ);
	free(err);
	free(phase);
	free(phase_dft);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe_dft(glwe_computed_dft);
	delete_module_info_p(module);
	delete_glwe_ct_params(params_glwe);
	delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Test glwe_secret_masking_dft. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking_dft, uniform_RnX_message)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	//! Parameters
	MODULE* module            = new_module_info_p(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKeyDFT* sk_dft      = new_glwe_secret_key_dft(NBASE, KBASE);
	GLWECiphertextDFT* glwe_computed_dft = new_glwe_dft(params_glwe);
	double* m_univ                = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* phase                = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyBivDFT* phase_dft         = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* m                    = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* err                  = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase_computed       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	double* phase_computed_univ   = calloc(NBASE, sizeof(double));

	//! Draws each input variable

	// Draws normaly in Rn[X] m_univ
	normal_random_vec(NBASE, m_univ, 1, NBASE, 0.0, 0.1);
	
	// Draws normaly the error
	normal_random_biv_poly(module, params_glwe, err);

	//! Computation with functions
	// Computes m_univ bivariate form
	univ_to_biv(params_glwe, m, m_univ);

	// Computes the final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking_dft(module, glwe_computed_dft, sk_dft, phase_dft);

	// Computes the computed phase in Rn[X]
	glwe_secret_demasking_dft(module, phase_computed, sk_dft, glwe_computed_dft);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff_1 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p];
		double diff_2 = m_univ[p] - round(m_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma + ldexp(1.0, -LBASE * KAPPABASE);

		// Using the triangle inequality, for each p, the difference should be smaller than :
		// |err_p| + |m_p - m_computed_p| Ie 3*sigma + 2^(-l*kappa)
		int cond = (diff_1 <= err_length && diff_1 >= -err_length) || (diff_2 <= err_length && diff_2 >= -err_length);

		if (!cond) big_error_count++;
	}

	// Asserts big_error_count <= 0.27*N
	cr_assert(big_error_count <= (int)(0.27 * NBASE), "The error should be greater than 3*sigma at most %ld times", (int)(0.27*NBASE));

	free(m);
	free(m_univ);
	free(err);
	free(phase);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe_dft(glwe_computed_dft);
	delete_module_info_p(module);
	delete_glwe_ct_params(params_glwe);
	delete_glwe_secret_key_dft(sk_dft);
}
