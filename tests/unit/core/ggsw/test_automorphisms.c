
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
#include "glwegadget_key.h"
#include "glwegadget_utils.h"
#include "maths_structures.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

PvdaParamTest(automorphism, generate_key, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEAutomorphismKey* auto_key  = new_automorphism_key(params_glwegadget);

	uniform_glwe_secret_key(module, sk, 3);
	glwe_sk_prepare(module, sk_prep, sk);

	int c = compute_automorphism_key(module, auto_key, sk_prep, 7);

	cr_assert(c == 0);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_automorphism_key(auto_key);
	DELETE_PVDA_PARAMS_GGSWGAD;
}

PvdaParamTest(automorphism, no_noise, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEAutomorphismKey* auto_key  = new_automorphism_key(params_glwegadget);

	PolyUnivTnX* m_univ_tnx     = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_observed_tnx = new_univ_tnx(params_glwe);
	PolyBiv* m_auto             = new_biv(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);
	GLWECiphertext* glwe_res    = new_glwe(params_glwe);
	GLWECiphertext* glwe_norm   = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 64);
	glwe_secret_encrypt_tnx(module, glwe_ct, sk_prep, m_univ_tnx);

	int auto_ps[] = {1, 3, 5, 7, 9, 13, 21, 29, 45, 91, -1, -3};

	for (int i = 0; i < sizeof(auto_ps) / sizeof(auto_ps[0]); ++i)
	{
		int auto_p = auto_ps[i];
		compute_automorphism_key(module, auto_key, sk_prep, auto_p);

		glwegadget_automorphism(module, glwe_res, auto_key, glwe_ct);
		normalize_glwe(module, glwe_norm, glwe_res);

		glwe_secret_decrypt(module, m_auto, sk_prep, glwe_norm);
		biv_to_univ_tnx(params_glwe, m_observed_tnx, m_auto);

		uint64_t nn = params_glwe->nn;
		pvda_znx_automorphism(module, auto_p, m_expected_tnx, m_univ_tnx);

		int64_t decomp_noise_bits = info_bits_half_prod(params_glwe, params_glwegadget);

		for (int p = 0; p < params_glwe->nn; ++p)
			assert_tnx_close_enough(m_observed_tnx[p], m_expected_tnx[p], decomp_noise_bits);
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_automorphism_key(auto_key);

	delete_univ_tnx(m_univ_tnx);
	delete_univ_tnx(m_observed_tnx);
	delete_univ_tnx(m_expected_tnx);
	delete_glwe(glwe_ct);
	delete_glwe(glwe_norm);
	delete_glwe(glwe_res);

	delete_biv(m_auto);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

PvdaParamTest(automorphism, noise, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEAutomorphismKey* auto_key  = new_automorphism_key(params_glwegadget);

	PolyUnivTnX* m_univ_tnx     = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_observed_tnx = new_univ_tnx(params_glwe);
	PolyBiv* m_auto             = new_biv(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);
	GLWECiphertext* glwe_res    = new_glwe(params_glwe);
	GLWECiphertext* glwe_norm   = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 64);
	glwe_secret_encrypt_tnx(module, glwe_ct, sk_prep, m_univ_tnx);

	int auto_ps[] = {1, 3, 5, 7, 9, 13, 21, 29, 45, 91, -1, -3};

	for (int i = 0; i < sizeof(auto_ps) / sizeof(auto_ps[0]); ++i)
	{
		int auto_p = auto_ps[i];
		compute_automorphism_key(module, auto_key, sk_prep, auto_p);

		glwegadget_automorphism(module, glwe_res, auto_key, glwe_ct);
		normalize_glwe(module, glwe_norm, glwe_res);

		glwe_secret_decrypt(module, m_auto, sk_prep, glwe_norm);
		biv_to_univ_tnx(params_glwe, m_observed_tnx, m_auto);

		//TODO: univ auto function
		uint64_t nn = params_glwe->nn;
		pvda_znx_automorphism(module, auto_p, m_expected_tnx, m_univ_tnx);

		// Not the real noise bound, but it seems like either the info_bits estimation
		// is more conservative than usually needed
		int64_t decomp_noise_bits = info_bits_half_prod(params_glwe, params_glwegadget);

		for (int p = 0; p < params_glwe->nn; ++p)
			assert_tnx_close_enough(m_observed_tnx[p], m_expected_tnx[p], decomp_noise_bits);
	}

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_automorphism_key(auto_key);

	delete_univ_tnx(m_univ_tnx);
	delete_univ_tnx(m_observed_tnx);
	delete_univ_tnx(m_expected_tnx);
	delete_glwe(glwe_ct);
	delete_glwe(glwe_norm);
	delete_glwe(glwe_res);

	delete_biv(m_auto);

	DELETE_PVDA_PARAMS_GGSWGAD;
}
