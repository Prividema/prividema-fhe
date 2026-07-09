#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw_arithmetic.h"
#include "core/glwe/glwe_arithmetic.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "ggsw_params.h"
#include "glwe_params.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"

/** The test is done without error, it is a proof of concept*/
PvdaParamTest(ggsw_external_product, without_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	//! Variance of the error's normal distributions
	params_glwe->fast_uniform_nb_bits = 0;
	sigma                             = 0;
	double err_length                 = glwe_bivariate_epsilon(params_glwe) + 3 * sigma + 3 * DBL_EPSILON;
	double critical_err_length        = glwe_bivariate_epsilon(params_glwe) + 5 * sigma + 5 * DBL_EPSILON;

	GLWESecretKey* sk_ggsw              = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_glwe_prep = alloc_glwe_secret_key_prepared(params_glwe);
	GGSWCiphertext* ggsw                = new_ggsw(params_ggsw);
	GLWECiphertext* glwe_tilde          = new_glwe(params_glwe);
	GLWECiphertext* ext_prod_observed   = new_glwe(params_glwe);
	PolyUniv* u_univ                    = new_univ(params_glwe);
	PolyBiv* m                          = new_biv(params_glwe);

	PolyBiv* phase_observed           = new_biv(params_glwe);
	PolyUnivRnX* um_observed_univ_rnx = new_univ_rnx(params_glwe);
	PolyUnivDFT* u_univ_dft           = new_univ_dft(module);
	PolyBivDFT* um_dft                = new_biv_dft(params_glwe);
	PolyBiv* um                       = new_biv(params_glwe);
	PolyUnivRnX* um_univ_rnx          = new_univ_rnx(params_glwe);

	uniform_glwe_secret_key(module, sk_ggsw, 3);
	glwe_sk_prepare(module, sk_glwe_prep, sk_ggsw);

	// Draws uniformly both messages
	uniform_pow2_random_pol_znx(module, u_univ, params_glwe->nn, params_glwe->kappa);
	uniform_random_biv_poly(module, params_glwe, m, 1);

	// Computation with function
	glwe_secret_encrypt_phase(module, glwe_tilde, sk_glwe_prep, m);
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
	ggsw_unprepared_external_product(module, ext_prod_observed, glwe_tilde, ggsw);
	normalize_glwe(module, ext_prod_observed, ext_prod_observed);
	glwe_secret_decrypt(module, phase_observed, sk_glwe_prep, ext_prod_observed);
	biv_to_univ_rnx(params_glwe, um_observed_univ_rnx, phase_observed);

	//Computes u*m manually
	univ_coefs_to_dft(module, u_univ_dft, u_univ);
	pvda_svp_apply_dft(module, um_dft, ggsw_params_l_tilde_a(params_ggsw), u_univ_dft, m);
	biv_dft_to_coefs(module, params_glwe, um, um_dft);
	pvda_vec_znx_normalize_base2k(module, params_glwe->kappa, um, um);
	biv_to_univ_rnx(params_glwe, um_univ_rnx, um);

	//! Asserts um_computed_univ(X) = u * m_univ
	pvda_assert_polynomial_distance(params_glwe, um_observed_univ_rnx, um_univ_rnx, err_length, critical_err_length);

	// Clean up
	delete_biv(m);
	delete_univ(u_univ);
	delete_univ_dft(module, u_univ_dft);
	delete_biv(phase_observed);
	delete_biv(um);
	delete_univ_rnx(um_univ_rnx);
	free(um_dft);
	delete_univ_rnx(um_observed_univ_rnx);

	delete_glwe(ext_prod_observed);
	delete_glwe(glwe_tilde);
	delete_ggsw(ggsw);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_prepared(sk_glwe_prep);

	DELETE_PVDA_PARAMS_GGSW;
}
