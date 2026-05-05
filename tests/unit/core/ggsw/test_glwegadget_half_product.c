#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>

#include "bivariate_polynomial.h"
#include "core/glwe/glwe_arithmetic.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

PvdaParamTest(glwegadget_half_product, without_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	sigma                             = 0;
	params_glwe->fast_uniform_nb_bits = 0;

	double biv_epsilon = glwe_bivariate_epsilon(params_glwe);
	double err_length =
	    params_glwe->nn * (2 * DBL_EPSILON + biv_epsilon) + 2 * glwe_params_l_a(params_glwe) * biv_epsilon;
	double critical_err_length =
	    params_glwe->nn * (3 * DBL_EPSILON + biv_epsilon) + 2 * glwe_params_l_a(params_glwe) * biv_epsilon;

	GLWESecretKey* sk                      = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft               = alloc_glwe_secret_key_dft(params_glwe);
	GLWECiphertext* glwe                   = new_glwe(params_glwe);
	GLWEGadgetCiphertext* glwegad          = new_glwegadget(params_glwegadget);
	GLWEGadgetCiphertextPrep* glwegad_prep = new_glwegadget_prep(params_glwegadget);
	PolyUniv* u_univ                       = new_univ(params_glwe);
	PolyUnivTnX* m_univ_tnx                = new_univ_tnx(params_glwe);
	PolyUnivRnX* m_univ_rnx                = new_univ_rnx(params_glwe);
	PolyUnivTnX* um_expected_tnx           = new_univ_tnx(params_glwe);
	PolyUnivRnX* um_expected_rnx           = new_univ_rnx(params_glwe);
	PolyBiv* um_observed                   = new_biv_poly(params_glwe);
	PolyUnivRnX* um_observed_rnx           = new_univ_rnx(params_glwe);
	PolyBiv* m                             = new_biv_poly(params_glwe);

	int nn    = params_glwe->nn;
	int k     = params_glwe->k;
	int kappa = params_glwe->kappa;

	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);
	uniform_random_pol_znx(u_univ, params_glwe->nn, 3);
	uniform_random_pol_znx(m_univ_tnx, params_glwe->nn, 62);
	univ_tnx_to_biv(params_glwe, m, m_univ_tnx, 0);

	memset(um_expected_tnx, 0, poly_univ_rnx_bytes(params_glwe));
	for (int i = 0; i < params_glwe->nn; ++i)
		for (int j = 0; j < params_glwe->nn; ++j)
		{
			if (i + j < params_glwe->nn)
				um_expected_tnx[(i + j) % params_glwe->nn] += (uint64_t)u_univ[i] * m_univ_tnx[j];
			else
				um_expected_tnx[(i + j) % params_glwe->nn] -= (uint64_t)u_univ[i] * m_univ_tnx[j];
		}

	univ_tnx_to_rnx(params_glwe, um_expected_rnx, um_expected_tnx);

	glwegadget_secret_encrypt(module, glwegad, sk_dft, u_univ);
	glwegadget_prepare(module, glwegad_prep, glwegad);
	glwegadget_half_prod(module, glwe, glwegad_prep, m);
	normalize_glwe(module, glwe, glwe);
	glwe_secret_decrypt(module, um_observed, sk_dft, glwe);
	biv_to_univ_rnx(params_glwe, um_observed_rnx, um_observed);

	pvda_assert_polynomial_distance(params_glwe, um_observed_rnx, um_expected_rnx, err_length, critical_err_length);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_glwe(glwe);
	delete_glwegadget(glwegad);
	delete_glwegadget_prep(glwegad_prep);
	delete_univ(u_univ);
	delete_univ_tnx(m_univ_tnx);
	delete_univ_rnx(m_univ_rnx);
	delete_univ_tnx(um_expected_tnx);
	delete_univ_rnx(um_expected_rnx);
	free(um_observed);
	delete_univ_rnx(um_observed_rnx);
	free(m);
	DELETE_PVDA_PARAMS_GGSWGAD;
}
