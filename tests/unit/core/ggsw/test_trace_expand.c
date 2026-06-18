

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "core/glwe/glwe_arithmetic.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "ggsw_utils.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

struct criterion_test_params trace_params_fn()
{
	static PvdaTstParams default_params[] = {
	    {.nn                        = 1024,
	     .k                         = 1,
	     .kappa                     = 4,
	     .ciphertext_nb_limbs       = 8l * 2,
	     .ciphertext_nb_limbs_tilde = 8l * 2,
	     .sigma                     = 0},  // toy params, let default sigma
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5,
	     .ciphertext_nb_limbs_tilde = 9l * 5,
	     .sigma                     = 0},  // k > 1 params
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5 - 1,
	     .ciphertext_nb_limbs_tilde = 9l * 5 - 1,
	     .sigma                     = 0},  // k > 1 l_a != l_b params
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 12l * 5,
	     .ciphertext_nb_limbs_tilde = 12l * 5,
	     .sigma                     = 0},  // k > 1 l_a != l_b params

	};

	return cr_make_param_array(PvdaTstParams, default_params, sizeof(default_params) / sizeof(default_params[0]));
}
PvdaParamTest(trace_expand, no_noise, trace_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUnivRnX* m_univ_rnx     = new_univ_rnx(params_glwe);
	PolyUnivRnX* tmp_rnx        = new_univ_rnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivRnX* m_observed_rnx = new_univ_rnx(params_glwe);
	PolyBiv* biv_tmp            = new_biv(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, params_glwe->nn, params_glwe->nn - 1};

	GLWEAutomorphismKSKCollection* ksks = new_automorphism_ksk_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKSK* ksk = new_automorphism_ksk(params_glwegadget);
		prepare_automorphism_key(module, ksk, sk_prep, p);
		glwegadget_ksk_collection_put_key(ksks, ksk, p);
	}
	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ_rnx, 0, poly_univ_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		rnx_random_vec(tmp_rnx, params_glwe);
		for (int i = 0; i < bund; ++i) m_univ_rnx[i] = tmp_rnx[i];

		glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

		GLWECiphertext** results = calloc(bund, sizeof(GLWECiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_glwe(params_glwe);
		}

		glwe_trace_expand(module, results, bund, glwe_ct, ksks);

		int a = 1;
		for (int i = 0; i < bund; ++i)
		{
			glwe_secret_decrypt(module, biv_tmp, sk_prep, results[i]);
			biv_to_univ_rnx(params_glwe, m_observed_rnx, biv_tmp);
			int64_t factor  = 1l << (next_pow2_log(bund));
			double expected = factor * m_univ_rnx[i];
			double actual   = m_observed_rnx[0];
			cr_assert(lt(dbl, rnx_torus_distance(expected, actual), 0.001));
			for (int p = 1; p < params_glwe->nn; ++p)
			{
				cr_assert(lt(dbl, rnx_torus_distance(0, m_observed_rnx[p]), 0.001));
			}
			delete_glwe(results[i]);
		}
		free(results);
	}

	delete_automorphism_ksk_collection(ksks, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_biv(biv_tmp);
	delete_univ_rnx(tmp_rnx);
	delete_univ_rnx(m_univ_rnx);
	delete_univ_rnx(m_observed_rnx);
	delete_univ_tnx(m_expected_tnx);
	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

struct criterion_test_params trace_params2_fn()
{
	static PvdaTstParams default_params[] = {
	    {.nn                        = 1024,
	     .k                         = 1,
	     .kappa                     = 4,
	     .ciphertext_nb_limbs       = 8l * 2,
	     .ciphertext_nb_limbs_tilde = 8l * 2,
	     .sigma                     = 0},  // toy params, let default sigma
	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5,
	     .ciphertext_nb_limbs_tilde = 9l * 5,
	     .sigma                     = 0},  // k > 1 params
	    {.nn                        = 256,
	     .k                         = 1,
	     .kappa                     = 19,
	     .ciphertext_nb_limbs       = 15l * 2,
	     .ciphertext_nb_limbs_tilde = 15l * 2,
	     .sigma                     = 0},  // k > 1 params

	    /*
	        {.nn                        = 1024,
	         .k                         = 4,
	         .kappa                     = 8,
	         .ciphertext_nb_limbs       = 9l * 5 - 1,
	         .ciphertext_nb_limbs_tilde = 9l * 5,
	         .sigma                     = 0},  // k > 1 l_a != l_b params
	     */

	};

	return cr_make_param_array(PvdaTstParams, default_params, sizeof(default_params) / sizeof(default_params[0]));
}

PvdaParamTest(ggsw_trace_expand, no_noise, trace_params2_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;
	double biv_epsilon                = glwe_bivariate_epsilon(params_glwe);
	double tst_epsilon                = DBL_EPSILON;
	double multiplier                 = params_glwe->nn;
	double max_err_length             = 2 * biv_epsilon;
	double critical_err_length        = 2 * biv_epsilon;

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUniv* m_univ        = new_univ(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, 1};

	GLWEAutomorphismKSKCollection* ksks = new_automorphism_ksk_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKSK* ksk = new_automorphism_ksk(params_glwegadget);
		prepare_automorphism_key(module, ksk, sk_prep, p);
		glwegadget_ksk_collection_put_key(ksks, ksk, p);
	}

	GGSWCiphertextPrep** ggsw_ksks = calloc(k, sizeof(GGSWCiphertext*));
	GGSWCiphertext* ggsw_tmp       = new_ggsw(params_ggsw);
	PolyUniv* neg_sk_i             = new_univ(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw_prep(params_ggsw);
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, ggsw_tmp, sk_prep, neg_sk_i);
		ggsw_prepare(module, ggsw_ksks[i], ggsw_tmp);
	}
	delete_univ(neg_sk_i);
	delete_ggsw(ggsw_tmp);

	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ, 0, poly_univ_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		uniform_random_pol_znx(m_univ, bund, 1);

		glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, bund);

		GGSWCiphertext** results = calloc(bund, sizeof(GGSWCiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_ggsw(params_ggsw);
		}

		packed_glwegadget_trace_expand_ggsw(module, results, bund, params_glwegadget->l_tilde, glwe_ct, ksks,
		                                    (const GGSWCiphertextPrep**)ggsw_ksks);

		PolyUniv* expected_b = new_univ(params_glwe);
		for (int b = 0; b < bund; ++b)
		{
			memset(expected_b, 0, poly_univ_bytes(params_glwe));
			expected_b[0] = m_univ[b];
			check_ggsw(module, results[b], sk_prep, expected_b, max_err_length, critical_err_length);
		}
	}

	delete_automorphism_ksk_collection(ksks, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

PvdaParamTest(glwegad2_trace_expand, no_noise, trace_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;
	double biv_epsilon                = glwe_bivariate_epsilon(params_glwe);
	double tst_epsilon                = DBL_EPSILON;
	double multiplier                 = params_glwe->nn;
	double max_err_length             = 3 * sigma + multiplier * tst_epsilon + 2 * biv_epsilon;
	double critical_err_length        = 5 * sigma + multiplier * tst_epsilon + 2 * biv_epsilon;

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);

	PolyUniv* m_univ        = new_univ(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	uniform_glwe_secret_key(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, 1};

	GLWEAutomorphismKSKCollection* ksks = new_automorphism_ksk_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKSK* ksk = new_automorphism_ksk(params_glwegadget);
		prepare_automorphism_key(module, ksk, sk_prep, p);
		glwegadget_ksk_collection_put_key(ksks, ksk, p);
	}

	GGSWCiphertext** ggsw_ksks = calloc(k, sizeof(GGSWCiphertext*));
	PolyUniv* neg_sk_i         = new_univ(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw(params_ggsw);
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_prepared_sk_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, ggsw_ksks[i], sk_prep, neg_sk_i);
	}
	delete_univ(neg_sk_i);

	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ, 0, poly_univ_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		uniform_random_pol_znx(m_univ, bund, 1);

		glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, bund);

		GLWEGadgetCiphertext** results = calloc(bund, sizeof(GGSWCiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_glwegadget(params_glwegadget);
		}

		packed_glwegadget_trace_expand(module, results, bund, params_glwegadget->l_tilde, glwe_ct, ksks);

		PolyUniv* expected_b = new_univ(params_glwe);
		for (int b = 0; b < bund; ++b)
		{
			memset(expected_b, 0, poly_univ_bytes(params_glwe));
			expected_b[0] = m_univ[b];
			check_glwegadget(module, results[b], sk_prep, expected_b, max_err_length, critical_err_length);
		}
	}

	delete_automorphism_ksk_collection(ksks, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}
