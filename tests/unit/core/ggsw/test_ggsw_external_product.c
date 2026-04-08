#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "rng.h"
#include "test_utils.h"
#include "univariate_polynomial.h"
#include "utils.h"

/** The test is done without error, it is a proof of concept*/
//TODO: is it, though? SIGMA_TILDEBASE is NOT 0
PvdaParamTest(ggsw_external_product, without_error, default_params_fn)
{
	INIT_PVDA_PARAMS_GGSW(param);

	//! Variance of the error's normal distributions
	params_glwe->sigma = 0;
	sigma              = 0;
	double err_length  = ldexp(1.0, -params_glwe->l * params_glwe->kappa) + 3 * sigma;

	GLWESecretKey* sk_ggsw            = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_glwe_dft     = alloc_glwe_secret_key_dft(params_glwe);
	GGSWCiphertext* ggsw              = new_ggsw(params_ggsw);
	GLWECiphertext* glwe_tilde        = new_glwe(params_glwe);
	GLWECiphertext* ext_prod_observed = new_glwe(params_glwe);
	PolyUniv* u_univ                  = new_univ(params_glwe);
	PolyBiv* m                        = new_biv_poly(params_glwe);

	PolyBiv* phase_observed           = new_biv_poly(params_glwe);
	PolyUnivRnX* um_observed_univ_RnX = new_univ_rnx(params_glwe);
	PolyUnivDFT* u_univ_dft           = new_univ_dft(module);
	PolyBivDFT* um_dft                = new_biv_poly_dft(params_glwe);
	PolyBiv* um                       = new_biv_poly(params_glwe);
	PolyUnivRnX* um_univ_RnX          = new_univ_rnx(params_glwe);

	// TODO: WHY?
	sk_ggsw->values[0] = 1;
	transform_glwe_secret_key_not_dft_to_dft(module, sk_glwe_dft, sk_ggsw);

	// Draws uniformly both messages
	uniform_random_pol_znx(u_univ, params_glwe->nn, params_glwe->kappa);
	uniform_random_biv_poly(params_glwe, m, 1);

	// Computation with function
	glwe_secret_encrypt_phase(module, glwe_tilde, sk_glwe_dft, m);
	ggsw_secret_encrypt(module, ggsw, sk_glwe_dft, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
	ggsw_external_product(module, ext_prod_observed, glwe_tilde, ggsw);
	normalize_glwe(module, ext_prod_observed, ext_prod_observed);
	glwe_secret_decrypt(module, phase_observed, sk_glwe_dft, ext_prod_observed);
	biv_to_univ_rnx(params_glwe, um_observed_univ_RnX, phase_observed);

	//Computes u*m manually
	univ_coefs_to_dft(module, u_univ_dft, u_univ);
	pvda_svp_apply_dft(module, um_dft, params_ggsw->l_tilde, u_univ_dft, m, params_ggsw->l_tilde, params_glwe->nn);
	univ_dft_to_coefs(module, um, um_dft);
	pvda_vec_znx_normalize_base2k(module, params_glwe->kappa, um, params_ggsw->l_tilde, params_glwe->nn, um,
	                              params_ggsw->l_tilde, params_glwe->nn);
	biv_to_univ_rnx(params_glwe, um_univ_RnX, um);

	//! Asserts um_computed_univ(X) = u * m_univ
	for (uint64_t p = 0; p < params_glwe->nn; p++)
	{
		double diff = torus_distance(um_univ_RnX[p], um_observed_univ_RnX[p]);
		int cond    = diff < err_length;

		cr_assert(cond, "Equality failed with um_computed_univ_RnX[%ld] = %lf and  um_univ_RnX[%ld] = %lf", p, p,
		          um_univ_RnX[p] - floor(um_univ_RnX[p]), p, um_observed_univ_RnX[p], err_length);
	}

	// Clean up
	free(m);
	delete_univ(u_univ);
	delete_univ_dft(u_univ_dft);
	free(phase_observed);
	free(um);
	delete_univ_rnx(um_univ_RnX);
	free(um_dft);
	delete_univ_rnx(um_observed_univ_RnX);

	delete_glwe(ext_prod_observed);
	delete_glwe(glwe_tilde);
	delete_ggsw(ggsw);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_dft(sk_glwe_dft);

	DELETE_PVDA_PARAMS_GGSW;
}
