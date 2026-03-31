#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "bivariate_polynomial.h"
#include "common/spqlios_alias.h"
#include "core/glwe/glwe.h"
#include "glwe_key.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

#define NBASE       1024
#define KBASE       1
#define KAPPABASE   4
#define NLIMBSBASE  (KBASE + 1) * 8
#define LBASE       NLIMBSBASE / (KBASE + 1)
#define SIGMABASE   -(LBASE / 2 + 1) * KAPPABASE

#define PROB_FACTOR 3

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

	// Since the message are drawn in Zn[X,Y], there is no decomposition error. Thus, the error should be smaller than 3*sigma 99.73% of the time
	double err_length = 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyBiv* m                           = new_biv_poly(params_glwe);
	PolyBiv* err                         = new_biv_poly(params_glwe);
	PolyUnivRnX* m_univ_RnX              = new_univ_rnx(params_glwe);
	PolyBiv* phase                       = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_computed        = new_glwe(params_glwe);
	PolyBiv* phase_computed              = new_biv_poly(params_glwe);
	PolyUnivRnX* phase_computed_univ_RnX = new_univ_rnx(params_glwe);

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// The input message uniformly drawn in Zn[X,Y]
	uniform_random_biv_poly(params_glwe, m, LBASE / 2);

	// The input error normaly drawn in Zn[X,Y]
	normal_random_biv_poly(params_glwe, err);

	//! Computations with functions
	// Computes the message in Tn[X] with the base-2Kappa
	biv_rnx_to_univ(params_glwe, m_univ_RnX, m);

	// The final phase = m + err
	add_biv_poly(module, params_glwe, phase, m, err);

	// Computes the bivGLWE ciphertext
	glwe_secret_encrypt_phase(module, glwe_computed, sk_dft, phase);

	// The computed phase in Rn[X]
	glwe_secret_decrypt(module, phase_computed, sk_dft, glwe_computed);

	// The computed phase in Rn[X]
	biv_rnx_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_computed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
	double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);
	//! Clean up
	free(phase_computed_univ_RnX);
	free(phase_computed);
	delete_glwe(glwe_computed);
	free(phase);
	free(m_univ_RnX);
	free(err);
	free(m);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
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

	// The message is drawn in Rn[X], there is a decomposition error of 2^{-kappa * l}. And the error should be smaller than 3*sigma 99.73% of the time.
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyBiv* m                           = new_biv_poly(params_glwe);
	PolyBiv* err                         = new_biv_poly(params_glwe);
	PolyUnivRnX* m_univ_RnX              = new_univ_rnx(params_glwe);
	PolyBiv* phase                       = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_computed        = new_glwe(params_glwe);
	PolyBiv* phase_computed              = new_biv_poly(params_glwe);
	PolyUnivRnX* phase_computed_univ_RnX = new_univ_rnx(params_glwe);
	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Draws normaly in Rn[X] m_univ_RnX
	normal_random_vec(m_univ_RnX, NBASE, 0.0, 0.1);

	// Draws normaly the error
	normal_random_biv_poly(params_glwe, err);

	//! Computation with functions
	// Computes m_univ_RnX bivariate form
	univ_rnx_to_biv(params_glwe, m, m_univ_RnX);

	// Computes the final phase = m + err
	add_biv_poly(module, params_glwe, phase, m, err);

	// Computes the bivGLWE ciphertext
	glwe_secret_encrypt_phase(module, glwe_computed, sk_dft, phase);

	// Computes the computed phase in Rn[X]
	glwe_secret_decrypt(module, phase_computed, sk_dft, glwe_computed);

	// The computed phase in Rn[X]
	biv_rnx_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_computed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
	double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);

	//! Clean up
	free(phase_computed_univ_RnX);
	free(phase_computed);
	delete_glwe(glwe_computed);
	free(phase);
	free(m);
	free(err);
	free(m_univ_RnX);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}
