

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "core/glwe/glwe_ciphertext.h"
#include "ggsw_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "ggsw_utils.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwegadget_arithmetic.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"
#include "glwegadget_utils.h"
#include "math_utils.h"
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

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);

	PolyUnivRnX* m_univ_rnx     = new_univ_rnx(params_glwe);
	PolyUnivRnX* tmp_rnx        = new_univ_rnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivRnX* m_observed_rnx = new_univ_rnx(params_glwe);
	PolyBiv* biv_tmp            = new_biv(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);

	glwe_sk_random_uniform(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, params_glwe->nn, params_glwe->nn - 1};

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* automorphism_key = new_automorphism_key(module, params_glwegadget);
		glwe_automorphism_key_generate(module, automorphism_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, automorphism_key, p);
	}
	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ_rnx, 0, univ_znx_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		rnx_random_vec(module, tmp_rnx, params_glwe);
		for (int i = 0; i < bund; ++i) m_univ_rnx[i] = tmp_rnx[i];

		glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

		GLWECiphertext** results = calloc(bund, sizeof(GLWECiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_glwe(params_glwe);
		}

		glwe_trace_expand(module, results, bund, glwe_ct, auto_key_collection);

		int a = 1;
		for (int i = 0; i < bund; ++i)
		{
			glwe_secret_decrypt(module, biv_tmp, sk_prep, results[i]);
			biv_to_univ_rnx(params_glwe, m_observed_rnx, biv_tmp, 0);
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

	delete_automorphism_key_collection(auto_key_collection, 1);

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

PvdaParamTest(hom_expand, no_noise, trace_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);

	PolyUnivRnX* m_univ_rnx     = new_univ_rnx(params_glwe);
	PolyUnivRnX* tmp_rnx        = new_univ_rnx(params_glwe);
	PolyUnivTnX* m_expected_tnx = new_univ_tnx(params_glwe);
	PolyUnivRnX* m_observed_rnx = new_univ_rnx(params_glwe);
	PolyBiv* biv_tmp            = new_biv(params_glwe);
	GLWECiphertext* glwe_ct     = new_glwe(params_glwe);

	glwe_sk_random_uniform(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, params_glwe->nn, params_glwe->nn - 1};

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* automorphism_key = new_automorphism_key(module, params_glwegadget);
		glwe_automorphism_key_generate(module, automorphism_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, automorphism_key, p);
	}
	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ_rnx, 0, univ_znx_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		rnx_random_vec(module, tmp_rnx, params_glwe);
		for (int i = 0; i < bund; ++i) m_univ_rnx[i] = tmp_rnx[i];

		glwe_secret_encrypt_rnx(module, glwe_ct, sk_prep, m_univ_rnx);

		GLWECiphertext** results = calloc(bund, sizeof(GLWECiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_glwe(params_glwe);
		}

		for (int i = 0; i < bund && i < 64; ++i)
			glwe_hom_trace(module, results[i], bund, i, glwe_ct, auto_key_collection);

		int a = 1;
		for (int i = 0; i < bund && i < 64; ++i)
		{
			glwe_secret_decrypt(module, biv_tmp, sk_prep, results[i]);
			biv_to_univ_rnx(params_glwe, m_observed_rnx, biv_tmp, 0);
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

	delete_automorphism_key_collection(auto_key_collection, 1);

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

	    {.nn                        = 1024,
	     .k                         = 4,
	     .kappa                     = 8,
	     .ciphertext_nb_limbs       = 9l * 5 - 1,
	     .ciphertext_nb_limbs_tilde = 9l * 5,
	     .sigma                     = 0},  // k > 1 l_a != l_b params

	};

	return cr_make_param_array(PvdaTstParams, default_params, sizeof(default_params) / sizeof(default_params[0]));
}

PvdaParamTest(ggsw_trace_expand, no_noise, trace_params2_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;
	double max_err_length             = 0.001;
	double critical_err_length        = 0.001;

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);

	PolyUnivZnx* m_univ        = new_univ_znx(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	glwe_sk_random_uniform(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, 1};

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* automorphism_key = new_automorphism_key(module, params_glwegadget);
		glwe_automorphism_key_generate(module, automorphism_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, automorphism_key, p);
	}

	GGSWCiphertextPrep** ggsw_ksks = calloc(k, sizeof(GGSWCiphertext*));
	GGSWCiphertext* ggsw_tmp       = new_ggsw(params_ggsw);
	PolyUnivZnx* neg_sk_i             = new_univ_znx(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw_prep(module, params_ggsw);
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_sk_prepared_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, ggsw_tmp, sk_prep, neg_sk_i);
		ggsw_prepare(module, ggsw_ksks[i], ggsw_tmp);
	}
	delete_univ_znx(neg_sk_i);
	delete_ggsw(ggsw_tmp);

	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ, 0, univ_znx_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		uniform_pow2_random_pol_znx(module, m_univ, bund, 1);

		glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, bund);

		GGSWCiphertext** results = calloc(bund, sizeof(GGSWCiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_ggsw(params_ggsw);
		}

		packed_glwegadget_trace_expand_ggsw(module, results, bund, params_glwegadget->l_tilde, glwe_ct,
		                                    auto_key_collection, (const GGSWCiphertextPrep**)ggsw_ksks);

		PolyUnivZnx* expected_b = new_univ_znx(params_glwe);
		for (int b = 0; b < bund; ++b)
		{
			memset(expected_b, 0, univ_znx_bytes(params_glwe));
			expected_b[0] = m_univ[b];
			check_ggsw(module, results[b], sk_prep, expected_b, max_err_length, critical_err_length);
		}
	}

	delete_automorphism_key_collection(auto_key_collection, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

PvdaParamTest(glwegad2_trace_expand, normal, trace_params_fn)
{
	INIT_PVDA_PARAMS_GGSWGAD(param);

	double biv_epsilon         = glwe_bivariate_epsilon(params_glwe);
	double max_err_length      = 0.001;
	double critical_err_length = 0.001;

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);

	PolyUnivZnx* m_univ        = new_univ_znx(params_glwe);
	GLWECiphertext* glwe_ct = new_glwe(params_glwe);

	glwe_sk_random_uniform(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);

	int bundled[] = {2, 4, 3, 5, 14, 1};

	GLWEAutomorphismKeyCollection* auto_key_collection = new_automorphism_key_collection(2 * params_glwe->nn);

	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* automorphism_key = new_automorphism_key(module, params_glwegadget);
		glwe_automorphism_key_generate(module, automorphism_key, sk_prep, p);
		glwegadget_key_collection_put_key(auto_key_collection, automorphism_key, p);
	}

	GGSWCiphertext** ggsw_ksks = calloc(k, sizeof(GGSWCiphertext*));
	PolyUnivZnx* neg_sk_i         = new_univ_znx(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw(params_ggsw);
		for (int p = 0; p < nn; ++p)
		{
			neg_sk_i[p] = -glwe_sk_prepared_extract_poly_coefs(sk_prep, i)[p];
		}
		ggsw_secret_encrypt(module, ggsw_ksks[i], sk_prep, neg_sk_i);
	}
	delete_univ_znx(neg_sk_i);

	for (int i = 0; i < sizeof(bundled) / sizeof(bundled[0]); ++i)
	{
		int bund = bundled[i];

		memset(m_univ, 0, univ_znx_bytes(params_glwe));

		// Get the message in univariate RnX form for expected result
		uniform_pow2_random_pol_znx(module, m_univ, bund, 1);

		glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, bund);

		GLWEGadgetCiphertext** results = calloc(bund, sizeof(GGSWCiphertext*));
		for (int i = 0; i < bund; ++i)
		{
			results[i] = new_glwegadget(params_glwegadget);
		}

		packed_glwegadget_trace_expand(module, results, bund, params_glwegadget->l_tilde, glwe_ct, auto_key_collection);

		PolyUnivZnx* expected_b = new_univ_znx(params_glwe);
		for (int b = 0; b < bund; ++b)
		{
			memset(expected_b, 0, univ_znx_bytes(params_glwe));
			expected_b[0] = m_univ[b];
			check_glwegadget(module, results[b], sk_prep, expected_b, max_err_length, critical_err_length);
		}
	}

	delete_automorphism_key_collection(auto_key_collection, 1);

	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);

	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

// =============================================
// |                                           |
// |     Error paths of prepared expansions    |
// |                                           |
// =============================================

// The prepared trace expansions free each unprepared element right after preparing it. If a later preparation
// fails, the cleanup block must not free the elements already deleted (double free).
// Failures are injected through the backend vtable: pvda_vmp_prepare_contiguous fails on its fail_at-th call.

#include "backend_private.h"

static PVDA_VMP_PREPARE_CONTIGUOUS_F* real_vmp_prepare_contiguous;
static int vmp_prepare_calls;
static int vmp_prepare_fail_at;  // <= 0: never fail

static int failing_vmp_prepare_contiguous(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat, const int64_t* mat,
                                          uint64_t nrows, uint64_t ncols)
{
	if (++vmp_prepare_calls == vmp_prepare_fail_at) return -1;
	return real_vmp_prepare_contiguous(module, pmat, mat, nrows, ncols);
}

static void install_failing_vmp_prepare(PvdaBackend* module)
{
	real_vmp_prepare_contiguous            = module->vt.pvda_vmp_prepare_contiguous;
	module->vt.pvda_vmp_prepare_contiguous = failing_vmp_prepare_contiguous;
}

static void restore_vmp_prepare(PvdaBackend* module) { module->vt.pvda_vmp_prepare_contiguous = real_vmp_prepare_contiguous; }

static PvdaTstParams error_path_params = {.nn                        = 256,
                                          .k                         = 1,
                                          .kappa                     = 19,
                                          .ciphertext_nb_limbs       = 15l * 2,
                                          .ciphertext_nb_limbs_tilde = 15l * 2,
                                          .sigma                     = 0};

#define ERROR_PATH_BUND 3

static GLWEAutomorphismKeyCollection* error_path_auto_keys(PvdaBackend* module, const GLWEParams* params_glwe,
                                                           const GLWEGadgetParams* params_glwegadget,
                                                           const GLWESecretKeyPrepared* sk_prep)
{
	GLWEAutomorphismKeyCollection* collection = new_automorphism_key_collection(2 * params_glwe->nn);
	for (uint64_t i = 1; (1ULL << i) <= params_glwe->nn; ++i)
	{
		int64_t p                     = (int64_t)params_glwe->nn / (1LL << (i - 1)) + 1;
		GLWEAutomorphismKey* auto_key = new_automorphism_key(module, params_glwegadget);
		cr_assert(eq(i32, glwe_automorphism_key_generate(module, auto_key, sk_prep, p), 0));
		glwegadget_key_collection_put_key(collection, auto_key, p);
	}
	return collection;
}

Test(glwegad_trace_expand_prepared, prepare_failure_no_double_free)
{
	PvdaTstParams* param = &error_path_params;
	INIT_PVDA_PARAMS_GGSWGAD(param);

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);
	PolyUnivZnx* m_univ               = new_univ_znx(params_glwe);
	GLWECiphertext* glwe_ct        = new_glwe(params_glwe);

	glwe_sk_random_uniform(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);
	GLWEAutomorphismKeyCollection* auto_keys = error_path_auto_keys(module, params_glwe, params_glwegadget, sk_prep);

	memset(m_univ, 0, univ_znx_bytes(params_glwe));
	uniform_pow2_random_pol_znx(module, m_univ, ERROR_PATH_BUND, 1);
	glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, ERROR_PATH_BUND);

	// results[0] must be pre-allocated (see the header), the others are allocated by the function
	GLWEGadgetCiphertextPrep* results[ERROR_PATH_BUND] = {0};
	results[0] = new_glwegadget_prep(module, params_glwegadget);

	install_failing_vmp_prepare(module);

	// Dry run: count the prepare calls, so the failure can target the second element of the preparation loop
	vmp_prepare_calls   = 0;
	vmp_prepare_fail_at = 0;
	cr_assert(eq(i32, packed_glwegadget_trace_expand_prepared(module, results, ERROR_PATH_BUND,
	                                                          params_glwegadget->l_tilde, glwe_ct, auto_keys),
	             0));
	int total_calls = vmp_prepare_calls;
	cr_assert(ge(i32, total_calls, ERROR_PATH_BUND));

	vmp_prepare_calls   = 0;
	vmp_prepare_fail_at = total_calls - ERROR_PATH_BUND + 2;
	cr_assert(eq(i32, packed_glwegadget_trace_expand_prepared(module, results, ERROR_PATH_BUND,
	                                                          params_glwegadget->l_tilde, glwe_ct, auto_keys),
	             -1));

	restore_vmp_prepare(module);

	for (int r = 0; r < ERROR_PATH_BUND; ++r) delete_glwegadget_prep(results[r]);
	delete_automorphism_key_collection(auto_keys, 1);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_univ_znx(m_univ);
	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}

Test(ggsw_trace_expand_prepared, prepare_failure_no_double_free)
{
	PvdaTstParams* param = &error_path_params;
	INIT_PVDA_PARAMS_GGSWGAD(param);

	uint64_t k  = params_glwe->k;
	uint64_t nn = params_glwe->nn;

	GLWESecretKey* sk              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = new_glwe_secret_key_prepared(module, params_glwe);
	PolyUnivZnx* m_univ               = new_univ_znx(params_glwe);
	GLWECiphertext* glwe_ct        = new_glwe(params_glwe);

	glwe_sk_random_uniform(module, sk, 1);
	glwe_sk_prepare(module, sk_prep, sk);
	GLWEAutomorphismKeyCollection* auto_keys = error_path_auto_keys(module, params_glwe, params_glwegadget, sk_prep);

	// GGSW encryptions of -sk_i, prepared
	GGSWCiphertextPrep** ggsw_ksks = calloc(k, sizeof(GGSWCiphertextPrep*));
	GGSWCiphertext* ggsw_tmp       = new_ggsw(params_ggsw);
	PolyUnivZnx* neg_sk_i             = new_univ_znx(params_glwe);
	for (uint64_t i = 0; i < k; ++i)
	{
		ggsw_ksks[i] = new_ggsw_prep(module, params_ggsw);
		for (uint64_t p = 0; p < nn; ++p) neg_sk_i[p] = -glwe_sk_prepared_extract_poly_coefs(sk_prep, i)[p];
		cr_assert(eq(i32, ggsw_secret_encrypt(module, ggsw_tmp, sk_prep, neg_sk_i), 0));
		cr_assert(eq(i32, ggsw_prepare(module, ggsw_ksks[i], ggsw_tmp), 0));
	}
	delete_univ_znx(neg_sk_i);
	delete_ggsw(ggsw_tmp);

	memset(m_univ, 0, univ_znx_bytes(params_glwe));
	uniform_pow2_random_pol_znx(module, m_univ, ERROR_PATH_BUND, 1);
	glwegadget_packed_secret_encrypt(module, glwe_ct, params_glwegadget, sk_prep, m_univ, ERROR_PATH_BUND);

	GGSWCiphertextPrep* results[ERROR_PATH_BUND] = {0};

	install_failing_vmp_prepare(module);

	// Dry run: count the prepare calls, so the failure can target the second element of the preparation loop
	vmp_prepare_calls   = 0;
	vmp_prepare_fail_at = 0;
	cr_assert(eq(i32, packed_glwegadget_trace_expand_ggsw_prepared(module, results, params_ggsw, ERROR_PATH_BUND,
	                                                               params_glwegadget->l_tilde, glwe_ct, auto_keys,
	                                                               (const GGSWCiphertextPrep**)ggsw_ksks),
	             0));
	int total_calls = vmp_prepare_calls;
	cr_assert(ge(i32, total_calls, ERROR_PATH_BUND));

	vmp_prepare_calls   = 0;
	vmp_prepare_fail_at = total_calls - ERROR_PATH_BUND + 2;
	cr_assert(eq(i32, packed_glwegadget_trace_expand_ggsw_prepared(module, results, params_ggsw, ERROR_PATH_BUND,
	                                                               params_glwegadget->l_tilde, glwe_ct, auto_keys,
	                                                               (const GGSWCiphertextPrep**)ggsw_ksks),
	             -1));

	restore_vmp_prepare(module);

	for (int r = 0; r < ERROR_PATH_BUND; ++r) delete_ggsw_prep(results[r]);
	for (uint64_t i = 0; i < k; ++i) delete_ggsw_prep(ggsw_ksks[i]);
	free(ggsw_ksks);
	delete_automorphism_key_collection(auto_keys, 1);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_univ_znx(m_univ);
	delete_glwe(glwe_ct);

	DELETE_PVDA_PARAMS_GGSWGAD;
}
