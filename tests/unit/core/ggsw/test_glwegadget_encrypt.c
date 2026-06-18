#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "ggsw_utils.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwe_transform_key.h"
#include "glwegadget_ciphertext.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

/**
 * @brief Tests ggsw_secret_encrpyt
 *
 */
PvdaParamTest(glwegadgetsenc, works, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	double biv_epsilon         = glwe_bivariate_epsilon(params_glwe);
	double err_length          = biv_epsilon + 3 * sigma + 3 * DBL_EPSILON;
	double critical_err_length = biv_epsilon + 5 * sigma + 5 * DBL_EPSILON;
	cr_log_info("error length = %e", err_length);

	GLWESecretKey* sk                = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep   = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEGadgetCiphertext* glwegadget = new_glwegadget(params_glwegadget);
	PolyUniv* m_univ                 = new_univ(params_glwe);

	// Draws the message
	uniform_glwe_secret_key(module, sk, 3);
	glwe_sk_prepare(module, sk_prep, sk);
	uniform_random_pol_znx(m_univ, params_glwe->nn, params_glwe->kappa);

	glwegadget_secret_encrypt(module, glwegadget, sk_prep, m_univ);

	check_glwegadget(module, glwegadget, sk_prep, m_univ, err_length, critical_err_length);

	// Clean up
	delete_glwe_secret_key(sk);
	delete_univ(m_univ);
	delete_glwegadget(glwegadget);
	delete_glwe_secret_key_prepared(sk_prep);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

PvdaParamTest(glwegadgetpacked_encrypt, works, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	double biv_epsilon         = glwe_bivariate_epsilon(params_glwe);
	double err_length          = biv_epsilon + 3 * sigma + 3 * DBL_EPSILON;
	double critical_err_length = biv_epsilon + 5 * sigma + 5 * DBL_EPSILON;
	cr_log_info("error length = %e", err_length);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	PolyUniv* m_univ               = new_univ(params_glwe);

	PolyBiv* phase_computed              = new_biv(params_glwe);
	PolyUnivRnX* phase_observed_univ_rnx = new_univ_rnx(params_glwe);
	PolyUnivRnX* phase_expected_univ_rnx = new_univ_rnx(params_glwe);

	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	const int D = 16;
	// Draws the message
	uniform_glwe_secret_key(module, sk, 3);
	glwe_sk_prepare(module, sk_prep, sk);
	memset(m_univ, 0, poly_univ_bytes(params_glwe));
	uniform_random_pol_znx(m_univ, D, params_glwe->kappa);

	glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, D);

	glwe_secret_decrypt(module, phase_computed, sk_prep, glwe_ct);
	biv_to_univ_rnx(params_glwe, phase_observed_univ_rnx, phase_computed);

	// Computes the expected result  m / 2^{kappa_tilde * i}
	memset(phase_expected_univ_rnx, 0, poly_univ_rnx_bytes(params_glwe));
	int64_t divlog = next_pow2_log(D * params_glwegadget->l_tilde);
	for (int i = 1; i <= params_glwegadget->l_tilde; ++i)
		for (uint64_t p = 0; p < D; p++)
			phase_expected_univ_rnx[(i - 1) * D + p] =
			    ldexp((double)m_univ[p], -(params_glwegadget->kappa_tilde * i + divlog));

	pvda_assert_polynomial_distance(params_glwe, phase_observed_univ_rnx, phase_expected_univ_rnx, err_length,
	                                critical_err_length);
	// Clean up
	delete_univ_rnx(phase_observed_univ_rnx);
	delete_biv(phase_computed);
	delete_univ_rnx(phase_expected_univ_rnx);
	delete_univ(m_univ);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GGSWGAD;
}
