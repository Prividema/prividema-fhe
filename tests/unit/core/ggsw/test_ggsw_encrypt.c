#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <math.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw_arithmetic.h"
#include "core/glwe/glwe_arithmetic.h"
#include "ggsw_ciphertext.h"
#include "ggsw_utils.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

/**
 * @brief Tests ggsw_secret_encrpyt
 *
 */
PvdaParamTest(ggsw_secret_encrypt, works, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	//Acceptable error:
	// 3/5 sigma from the Normal error we add for encryption
	// N * DBL_EPSILON from N-len polynomial multiplication
	// 2 * biv_epsilon for bivariate to/from conversion
	double biv_epsilon         = glwe_bivariate_epsilon(params_glwe);
	double tst_epsilon         = DBL_EPSILON;
	double multiplier          = params_glwe->nn;
	double max_err_length      = 3 * sigma + tst_epsilon + 2 * biv_epsilon;
	double critical_err_length = 5 * sigma + tst_epsilon + 2 * biv_epsilon;

	GLWESecretKey* sk              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GGSWCiphertext* ggsw           = new_ggsw(params_ggsw);
	PolyUniv* m_univ               = new_univ(params_glwe);
	// Variables to compute the phase of each ggsw's row

	// Draws uniformly in (Zn[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	glwe_sk_prepare(module, sk_prep, sk);
	uniform_pow2_random_pol_znx(m_univ, params_glwe->nn, params_glwe->kappa);

	ggsw_secret_encrypt(module, ggsw, sk_prep, m_univ);

	check_ggsw(module, ggsw, sk_prep, m_univ, max_err_length, critical_err_length);

	// Clean up
	delete_univ(m_univ);
	delete_ggsw(ggsw);
	delete_glwe_secret_key_prepared(sk_prep);
	delete_glwe_secret_key(sk);

	DELETE_PVDA_PARAMS_GGSW;
}
