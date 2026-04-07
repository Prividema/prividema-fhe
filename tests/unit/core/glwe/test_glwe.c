#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "bivariate_polynomial.h"
#include "common/spqlios_alias.h"
#include "core/glwe/glwe.h"
#include "glwe_key.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"
#include "ututils.h"

PvdaTstParams params = {1024, 1, 4, 8, 0, 0};

#define PROB_FACTOR 3

//! GGWS PART (begin)

/**
 * @brief Test glwe_secret_masking. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking, small_error)
{
	INIT_PVDA_PARAMS_GLWE(&params);
	double err_length = 3 * sigma;

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(params_glwe);
	PolyBiv* m                           = new_biv_poly(params_glwe);
	PolyBiv* err                         = new_biv_poly(params_glwe);
	PolyUnivRnX* m_univ_RnX              = new_univ_rnx(params_glwe);
	PolyBiv* phase                       = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_observed        = new_glwe(params_glwe);
	PolyBiv* phase_observed              = new_biv_poly(params_glwe);
	PolyUnivRnX* phase_observed_univ_RnX = new_univ_rnx(params_glwe);

	//Draw key and message
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);
	uniform_random_biv_poly(params_glwe, m, params_glwe->l / 2);

	// Get the message in univariate RnX form for expected result
	biv_to_univ_rnx(params_glwe, m_univ_RnX, m);

	// Encrypt and decrypt the message
	normal_random_biv_poly(params_glwe, err);
	add_biv_poly(module, params_glwe, phase, m, err);
	glwe_secret_encrypt_phase(module, glwe_observed, sk_dft, phase);
	glwe_secret_decrypt(module, phase_observed, sk_dft, glwe_observed);
	biv_to_univ_rnx(params_glwe, phase_observed_univ_RnX, phase_observed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_observed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * params_glwe->nn);
	double proba  = binomial_tail(params_glwe->nn, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);
	//! Clean up
	free(phase_observed_univ_RnX);
	free(phase_observed);
	delete_glwe(glwe_observed);
	free(phase);
	free(m_univ_RnX);
	free(err);
	free(m);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);

	DELETE_PVDA_PARAMS_GLWE;
}

/**
 * @brief Test glwe_encrypt_priv. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking, uniform_RnX_message)
{
	INIT_PVDA_PARAMS_GLWE(&params);

	double err_length = ldexp(1.0, -params_glwe->l * params_glwe->kappa) + 3 * sigma;

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(params_glwe);
	PolyBiv* m                           = new_biv_poly(params_glwe);
	PolyBiv* err                         = new_biv_poly(params_glwe);
	PolyUnivRnX* m_univ_RnX              = new_univ_rnx(params_glwe);
	PolyBiv* phase                       = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_observed        = new_glwe(params_glwe);
	PolyBiv* phase_observed              = new_biv_poly(params_glwe);
	PolyUnivRnX* phase_observed_univ_RnX = new_univ_rnx(params_glwe);

	//Draw message (in RnX) and key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);
	normal_random_vec(m_univ_RnX, params_glwe->nn, 0.0, 0.1);

	//Encrypt said message
	univ_rnx_to_biv(params_glwe, m, m_univ_RnX);
	normal_random_biv_poly(params_glwe, err);
	add_biv_poly(module, params_glwe, phase, m, err);
	glwe_secret_encrypt_phase(module, glwe_observed, sk_dft, phase);
	glwe_secret_decrypt(module, phase_observed, sk_dft, glwe_observed);
	biv_to_univ_rnx(params_glwe, phase_observed_univ_RnX, phase_observed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_observed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * params_glwe->nn);
	double proba  = binomial_tail(params_glwe->nn, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);

	//! Clean up
	free(phase_observed_univ_RnX);
	free(phase_observed);
	delete_glwe(glwe_observed);
	free(phase);
	free(m);
	free(err);
	free(m_univ_RnX);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);

	DELETE_PVDA_PARAMS_GLWE;
}
