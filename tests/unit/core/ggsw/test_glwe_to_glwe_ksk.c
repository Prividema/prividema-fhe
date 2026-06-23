
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

PvdaParamTest(glwe_ksk, no_noise, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk                  = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep     = alloc_glwe_secret_key_prepared(params_glwe);
	GLWESecretKey* new_sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* new_sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEAutomorphismKSK* auto_ksk      = new_automorphism_ksk(params_glwegadget);

	PolyUnivTnX* m_univ_tnx     = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_observed_tnx = new_univ_tnx(params_glwe);
	PolyBiv* m_auto             = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);
	GLWECiphertext* glwe_res    = new_glwe(params_glwe);
	GLWECiphertext* glwe_norm   = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_glwe_secret_key(module, new_sk, 1);
	glwe_sk_prepare(module, new_sk_prep, new_sk);

	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 64);
	glwe_secret_encrypt_tnx(module, glwe_ct, sk_prep, m_univ_tnx);

	prepare_ksk(module, auto_ksk, new_sk_prep, sk_prep);

	glwe_to_glwe_keyswitch(module, glwe_res, auto_ksk, glwe_ct);
	normalize_glwe(module, glwe_norm, glwe_res);

	glwe_secret_decrypt(module, m_auto, new_sk_prep, glwe_norm);
	biv_to_univ_tnx(params_glwe, m_observed_tnx, m_auto);

	uint64_t nn = params_glwe->nn;

	int64_t decomp_noise_bits = info_bits_half_prod(params_glwe, params_glwegadget);

	for (int p = 0; p < params_glwe->nn; ++p)
		assert_tnx_close_enough(m_observed_tnx[p], m_univ_tnx[p], decomp_noise_bits);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_automorphism_ksk(auto_ksk);

	delete_univ_tnx(m_univ_tnx);
	delete_univ_tnx(m_observed_tnx);
	delete_glwe(glwe_ct);
	delete_glwe(glwe_norm);
	delete_glwe(glwe_res);

	delete_biv(m_auto);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

PvdaParamTest(glwe_ksk, noise, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	GLWESecretKey* sk                  = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep     = alloc_glwe_secret_key_prepared(params_glwe);
	GLWESecretKey* new_sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* new_sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GLWEAutomorphismKSK* auto_ksk      = new_automorphism_ksk(params_glwegadget);

	PolyUnivTnX* m_univ_tnx     = new_univ_tnx(params_glwe);
	PolyUnivTnX* m_observed_tnx = new_univ_tnx(params_glwe);
	PolyBiv* m_auto             = new_biv_poly(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);
	GLWECiphertext* glwe_res    = new_glwe(params_glwe);
	GLWECiphertext* glwe_norm   = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	uniform_glwe_secret_key(module, new_sk, 1);
	glwe_sk_prepare(module, new_sk_prep, new_sk);

	uniform_random_pol_znx((PolyUniv*)m_univ_tnx, params_glwe->nn, 64);
	glwe_secret_encrypt_tnx(module, glwe_ct, sk_prep, m_univ_tnx);

	prepare_ksk(module, auto_ksk, new_sk_prep, sk_prep);

	glwe_to_glwe_keyswitch(module, glwe_res, auto_ksk, glwe_ct);
	normalize_glwe(module, glwe_norm, glwe_res);

	glwe_secret_decrypt(module, m_auto, new_sk_prep, glwe_norm);
	biv_to_univ_tnx(params_glwe, m_observed_tnx, m_auto);

	uint64_t nn = params_glwe->nn;

	//Not the real bound (it should be looser), but it is enough for now
	int64_t decomp_noise_bits = info_bits_half_prod(params_glwe, params_glwegadget);

	for (int p = 0; p < params_glwe->nn; ++p)
		assert_tnx_close_enough(m_observed_tnx[p], m_univ_tnx[p], decomp_noise_bits);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_automorphism_ksk(auto_ksk);

	delete_univ_tnx(m_univ_tnx);
	delete_univ_tnx(m_observed_tnx);
	delete_glwe(glwe_ct);
	delete_glwe(glwe_norm);
	delete_glwe(glwe_res);

	delete_biv(m_auto);

	DELETE_PVDA_PARAMS_GGSWGAD;
}
