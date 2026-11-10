#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <float.h>
#include <math.h>

#include "backend_arithmetic.h"
#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw_arithmetic.h"
#include "core/glwe/glwe_arithmetic.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_key.h"
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
	double eps                        = glwe_bivariate_epsilon(params_glwe);

	// Even with sigma = 0, truncating the GGSW bodies to l_b limbs leaves a rounding noise of
	// variance eps^2/12 in every row. The external product sums it over ggsw_num_rows * nn terms,
	// each multiplied by a decomposition digit of variance 4^kappa/12.
	double sigma_ext = eps * sqrt((double)ggsw_num_rows(params_ggsw) * (double)params_glwe->nn *
	                              (ldexp(1.0, 2 * (int)params_glwe->kappa) / 12) / 12);

	double err_length          = eps + 3 * sigma_ext + 3 * DBL_EPSILON;
	double critical_err_length = eps + 5 * sigma_ext + 5 * DBL_EPSILON;

	GLWESecretKey* sk_ggsw              = new_glwe_secret_key(params_glwe);
	GLWESecretKeyPrepared* sk_glwe_prep = new_glwe_secret_key_prepared(module, params_glwe);
	GGSWCiphertext* ggsw                = new_ggsw(params_ggsw);
	GLWECiphertext* glwe_tilde          = new_glwe(params_glwe);
	GLWECiphertext* ext_prod_observed   = new_glwe(params_glwe);
	PolyUnivZnx* u_univ                    = new_univ_znx(params_glwe);
	PolyBiv* m                          = new_biv(params_glwe);

	PolyBiv* phase_observed           = new_biv(params_glwe);
	PolyUnivRnX* um_observed_univ_rnx = new_univ_rnx(params_glwe);
	BE_PolyUnivSvpPrep* u_univ_prep       = new_univ_svp_prep(module);
	BE_PolyBivDFT* um_dft             = new_biv_dft(module, params_glwe);
	PolyBiv* um                       = new_biv(params_glwe);
	PolyUnivRnX* um_univ_rnx          = new_univ_rnx(params_glwe);

	glwe_sk_random_uniform(module, sk_ggsw, 3);
	glwe_sk_prepare(module, sk_glwe_prep, sk_ggsw);

	// Draws uniformly both messages
	uniform_pow2_random_pol_znx(module, u_univ, params_glwe->nn, params_glwe->kappa);
	biv_random_uniform(module, params_glwe, m, 1);

	// Computation with function
	glwe_secret_encrypt_phase(module, glwe_tilde, sk_glwe_prep, m);
	ggsw_secret_encrypt(module, ggsw, sk_glwe_prep, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
	ggsw_external_product_unprepared(module, ext_prod_observed, glwe_tilde, ggsw);
	glwe_normalize(module, ext_prod_observed, ext_prod_observed);
	glwe_secret_decrypt(module, phase_observed, sk_glwe_prep, ext_prod_observed);
	biv_to_univ_rnx(params_glwe, um_observed_univ_rnx, phase_observed, 0);

	//Computes u*m manually
	univ_znx_to_svp_prep(module, u_univ_prep, u_univ);
	pvda_svp_apply_dft(module, um_dft, ggsw_params_l_tilde_a(params_ggsw), u_univ_prep, m);
	biv_from_dft(module, params_glwe, um, um_dft);
	biv_normalize(module, params_glwe->kappa, um, um);
	biv_to_univ_rnx(params_glwe, um_univ_rnx, um, 0);

	//! Asserts um_computed_univ(X) = u * m_univ
	pvda_assert_polynomial_distance(params_glwe, um_observed_univ_rnx, um_univ_rnx, err_length, critical_err_length);

	// Clean up
	delete_biv(m);
	delete_univ_znx(u_univ);
	delete_univ_svp_prep(module, u_univ_prep);
	delete_biv(phase_observed);
	delete_biv(um);
	delete_univ_rnx(um_univ_rnx);
	delete_biv_dft(module, um_dft);
	delete_univ_rnx(um_observed_univ_rnx);

	delete_glwe(ext_prod_observed);
	delete_glwe(glwe_tilde);
	delete_ggsw(ggsw);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_prepared(sk_glwe_prep);

	DELETE_PVDA_PARAMS_GGSW;
}
