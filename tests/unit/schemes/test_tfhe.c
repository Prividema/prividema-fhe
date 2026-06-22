
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw_arithmetic.h"
#include "core/glwe/glwe_arithmetic.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "rng.h"
#include "test_utils.h"
#include "tfhe.h"
#include "univariate_polynomial.h"

/** The test is done without error, it is a proof of concept*/
PvdaParamTest(tfhe_cmux_unprepared, without_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	//! Variance of the error's normal distributions
	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk_ggsw              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_glwe_prep = alloc_glwe_secret_key_prepared(params_glwe);

	GGSWCiphertext* ggsw = new_ggsw(params_ggsw);

	PolyUnivTnX* m1 = new_univ_tnx(params_glwe);
	PolyUnivTnX* m2 = new_univ_tnx(params_glwe);

	GLWECiphertext* glwe1 = new_glwe(params_glwe);
	GLWECiphertext* glwe2 = new_glwe(params_glwe);
	GLWECiphertext* res   = new_glwe(params_glwe);

	PolyBiv* res_biv     = new_biv(params_glwe);
	PolyUnivTnX* res_tnx = new_univ_tnx(params_glwe);

	PolyUniv* m_sel = new_univ(params_glwe);

	uniform_glwe_secret_key(module, sk_ggsw, 3);
	glwe_sk_prepare(module, sk_glwe_prep, sk_ggsw);

	uniform_random_pol_znx(m1, params_glwe->nn, 62);
	uniform_random_pol_znx(m2, params_glwe->nn, 62);

	glwe_secret_encrypt_tnx(module, glwe1, sk_glwe_prep, m1);
	glwe_secret_encrypt_tnx(module, glwe2, sk_glwe_prep, m2);

	//Case 0:
	memset(m_sel, 0, poly_univ_bytes(params_glwe));
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, m_sel);

	tfhe_cmux_unprepared(module, res, glwe1, glwe2, ggsw, 1);

	glwe_secret_decrypt(module, res_biv, sk_glwe_prep, res);

	biv_to_univ_tnx(params_glwe, res_tnx, res_biv);

	int decomp_noise_bits = 15;  //TODO: put a real threshold
	for (int p = 0; p < params_glwe->nn; ++p) assert_tnx_close_enough(res_tnx[p], m1[p], decomp_noise_bits);

	//Case 1:
	//
	//
	m_sel[0] = 1;
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, m_sel);

	tfhe_cmux_unprepared(module, res, glwe1, glwe2, ggsw, 1);

	glwe_secret_decrypt(module, res_biv, sk_glwe_prep, res);

	biv_to_univ_tnx(params_glwe, res_tnx, res_biv);

	for (int p = 0; p < params_glwe->nn; ++p) assert_tnx_close_enough(res_tnx[p], m2[p], decomp_noise_bits);

	delete_univ_tnx(m1);
	delete_univ_tnx(m2);

	delete_glwe(glwe1);
	delete_glwe(glwe2);
	delete_glwe(res);

	delete_ggsw(ggsw);

	delete_biv(res_biv);
	delete_univ_tnx(res_tnx);
	delete_univ(m_sel);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_prepared(sk_glwe_prep);

	DELETE_PVDA_PARAMS_GGSW;
}

/** The test is done without error, it is a proof of concept*/
PvdaParamTest(tfhe_cmux_prepared, without_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	//! Variance of the error's normal distributions
	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk_ggsw              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_glwe_prep = alloc_glwe_secret_key_prepared(params_glwe);

	GGSWCiphertext* ggsw          = new_ggsw(params_ggsw);
	GGSWCiphertextPrep* ggsw_prep = new_ggsw_prep(params_ggsw);

	PolyUnivTnX* m1 = new_univ_tnx(params_glwe);
	PolyUnivTnX* m2 = new_univ_tnx(params_glwe);

	GLWECiphertext* glwe1 = new_glwe(params_glwe);
	GLWECiphertext* glwe2 = new_glwe(params_glwe);
	GLWECiphertext* res   = new_glwe(params_glwe);

	PolyBiv* res_biv     = new_biv(params_glwe);
	PolyUnivTnX* res_tnx = new_univ_tnx(params_glwe);

	PolyUniv* m_sel = new_univ(params_glwe);

	uniform_glwe_secret_key(module, sk_ggsw, 3);
	glwe_sk_prepare(module, sk_glwe_prep, sk_ggsw);

	uniform_random_pol_znx(m1, params_glwe->nn, 62);
	uniform_random_pol_znx(m2, params_glwe->nn, 62);

	glwe_secret_encrypt_tnx(module, glwe1, sk_glwe_prep, m1);
	glwe_secret_encrypt_tnx(module, glwe2, sk_glwe_prep, m2);

	//Case 0:
	memset(m_sel, 0, poly_univ_bytes(params_glwe));
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, m_sel);

	ggsw_prepare(module, ggsw_prep, ggsw);
	tfhe_cmux(module, res, glwe1, glwe2, ggsw_prep, 1);

	glwe_secret_decrypt(module, res_biv, sk_glwe_prep, res);

	biv_to_univ_tnx(params_glwe, res_tnx, res_biv);

	int decomp_noise_bits = 15;  //TODO: put a real threshold
	for (int p = 0; p < params_glwe->nn; ++p) assert_tnx_close_enough(res_tnx[p], m1[p], decomp_noise_bits);

	//Case 1:
	//
	//
	m_sel[0] = 1;
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, m_sel);

	ggsw_prepare(module, ggsw_prep, ggsw);

	tfhe_cmux(module, res, glwe1, glwe2, ggsw_prep, 1);

	glwe_secret_decrypt(module, res_biv, sk_glwe_prep, res);

	biv_to_univ_tnx(params_glwe, res_tnx, res_biv);

	for (int p = 0; p < params_glwe->nn; ++p) assert_tnx_close_enough(res_tnx[p], m2[p], decomp_noise_bits);

	delete_univ_tnx(m1);
	delete_univ_tnx(m2);

	delete_glwe(glwe1);
	delete_glwe(glwe2);
	delete_glwe(res);

	delete_ggsw(ggsw);
	delete_ggsw_prep(ggsw_prep);

	delete_biv(res_biv);
	delete_univ_tnx(res_tnx);
	delete_univ(m_sel);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_prepared(sk_glwe_prep);

	DELETE_PVDA_PARAMS_GGSW;
}

/** The test is done without error, it is a proof of concept*/
GLWESecretKeyPrepared* dbg2k = NULL;

PvdaParamTest(tfhe_cmux_tree, without_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	//! Variance of the error's normal distributions
	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;

	GLWESecretKey* sk_ggsw              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_glwe_prep = alloc_glwe_secret_key_prepared(params_glwe);
	dbg2k                               = sk_glwe_prep;

	GGSWCiphertext* ggsw           = new_ggsw(params_ggsw);
	GGSWCiphertextPrep* ggsw_prep0 = new_ggsw_prep(params_ggsw);
	GGSWCiphertextPrep* ggsw_prep1 = new_ggsw_prep(params_ggsw);

	int num_msgs     = 64;
	int selections[] = {0, 1, 2, 3, 5, 7, 8, 13, 17, 31, 32, 33, 34};
	int msgs_log2    = 6;

	PolyUnivTnX* ms[num_msgs];

	GLWECiphertext* glwes[num_msgs];
	GLWECiphertext* res = new_glwe(params_glwe);

	PolyBiv* res_biv     = new_biv(params_glwe);
	PolyUnivTnX* res_tnx = new_univ_tnx(params_glwe);

	PolyUniv* m_sel = new_univ(params_glwe);

	uniform_glwe_secret_key(module, sk_ggsw, 2);
	glwe_sk_prepare(module, sk_glwe_prep, sk_ggsw);

	for (int i = 0; i < num_msgs; ++i)
	{
		ms[i]    = new_univ_tnx(params_glwe);
		glwes[i] = new_glwe(params_glwe);
		//uniform_random_pol_znx(ms[i], params_glwe->nn, 62);
		ms[i][0] = (uint64_t)i << 57;
		ms[i][2] = (uint64_t)(1 + i) << 32;
		ms[i][1] = (uint64_t)(13 + i) << 32;
		glwe_secret_encrypt_tnx(module, glwes[i], sk_glwe_prep, ms[i]);
	}

	// Prepare GGSW 0 and 1
	memset(m_sel, 0, poly_univ_bytes(params_glwe));
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, m_sel);

	ggsw_prepare(module, ggsw_prep0, ggsw);

	m_sel[0] = 1;
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, m_sel);

	ggsw_prepare(module, ggsw_prep1, ggsw);

	GGSWCiphertextPrep* selector[msgs_log2];
	for (int i = 0; i < sizeof(selections) / sizeof(selections[0]); ++i)
	{
		int sel = selections[i];
		for (int j = 0; j < msgs_log2; ++j)
		{
			uint64_t v = sel % 2;
			if (v)
			{
				selector[j] = ggsw_prep1;
			}
			else
			{
				selector[j] = ggsw_prep0;
			}
			v >>= 1;
		}

		tfhe_cmux_tree(module, res, (const GLWECiphertext**)glwes, num_msgs, (const GGSWCiphertextPrep**)selector,
		               msgs_log2, 0);
		glwe_secret_decrypt(module, res_biv, sk_glwe_prep, res);
		biv_to_univ_tnx(params_glwe, res_tnx, res_biv);

		int decomp_noise_bits = 5;  //TODO: put a real threshold
		for (int p = 0; p < params_glwe->nn; ++p) assert_tnx_close_enough(res_tnx[p], ms[sel][p], decomp_noise_bits);
	}

	delete_glwe(res);

	delete_ggsw(ggsw);
	delete_ggsw_prep(ggsw_prep0);
	delete_ggsw_prep(ggsw_prep1);

	delete_biv(res_biv);
	delete_univ_tnx(res_tnx);
	delete_univ(m_sel);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_prepared(sk_glwe_prep);

	DELETE_PVDA_PARAMS_GGSW;
}
