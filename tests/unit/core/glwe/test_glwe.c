#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

//! GGWS PART (begin)

/**
 * @brief Test glwe_secret_masking. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
PvdaParamTest(glwe_secret_masking, small_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	double biv_epsilon         = glwe_bivariate_epsilon(params_glwe);
	double err_length          = 3 * sigma + 2 * (biv_epsilon + DBL_EPSILON);
	double critical_err_length = 5 * sigma + 2 * (biv_epsilon + DBL_EPSILON);

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
	uniform_random_biv_poly(params_glwe, m, glwe_params_l_a(params_glwe));

	// Get the message in univariate RnX form for expected result
	biv_to_univ_rnx(params_glwe, m_univ_RnX, m);

	// Encrypt and decrypt the message
	normal_random_biv_poly(params_glwe, err);
	add_biv_poly(module, params_glwe, phase, m, err);
	glwe_secret_encrypt_phase(module, glwe_observed, sk_dft, phase);
	glwe_secret_decrypt(module, phase_observed, sk_dft, glwe_observed);
	biv_to_univ_rnx(params_glwe, phase_observed_univ_RnX, phase_observed);

	pvda_assert_polynomial_distance(params_glwe, phase_observed_univ_RnX, m_univ_RnX, err_length, critical_err_length);

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
PvdaParamTest(glwe_secret_masking, uniform_RnX_message, default_params_fn)
{
	INIT_PVDA_PARAMS_GLWE(param);

	double biv_epsilon         = glwe_bivariate_epsilon(params_glwe);
	double err_length          = 3 * sigma + 2 * (biv_epsilon + DBL_EPSILON);
	double critical_err_length = 5 * sigma + 2 * (biv_epsilon + DBL_EPSILON);

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
	univ_rnx_to_biv(params_glwe, m, m_univ_RnX, 0);
	normal_random_biv_poly(params_glwe, err);
	add_biv_poly(module, params_glwe, phase, m, err);
	glwe_secret_encrypt_phase(module, glwe_observed, sk_dft, phase);
	glwe_secret_decrypt(module, phase_observed, sk_dft, glwe_observed);
	biv_to_univ_rnx(params_glwe, phase_observed_univ_RnX, phase_observed);

	pvda_assert_polynomial_distance(params_glwe, phase_observed_univ_RnX, m_univ_RnX, err_length, critical_err_length);

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
