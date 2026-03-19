#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
#include "ggsw_params.h"
#include "glwegad.h"
#include "glwegad_ciphertext.h"
#include "rng.h"
#include "utils.h"

#define NBASE            8
#define KBASE            1
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 2
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        -(LBASE / 2 + 1) * KAPPABASE

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 2
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)
#define SIGMA_TILDEBASE  -3

Test(ggsw_external_product, without_error)
{
	//! Variance of the error's normal distributions
	double sigma       = 0;
	double sigma_tilde = 0;
	double err_length  = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;

	//! Parameters
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GLWEParams* params_glwe_tilde =
	    new_glwe_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
	GLWEGadParams* params_glwegad = new_glwegad_params(params_glwe, KAPPA_TILDEBASE, L_TILDEBASE);
	MODULE* module                = new_module_info(NBASE, FFT64);

	//! Variables
	GLWESecretKey* sk_ggsw            = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_glwe_dft     = alloc_glwe_secret_key_dft(NBASE, KBASE);
	GLWEGadCiphertext* glwegad_ct     = new_glwegad(params_glwegad);
	GLWECiphertext* glwe_tilde        = new_glwe(params_glwe_tilde);
	GLWECiphertext* ext_prod_computed = new_glwe(params_glwe);
	PolyUniv* u_univ                  = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* m                        = malloc(poly_biv_bytes(params_glwe_tilde));

	//! Variables to compute the result phase of the external product
	PolyBiv* phase_computed           = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* um_computed_univ_RnX = calloc(NBASE, sizeof(double));
	PolyUnivDFT* u_univ_dft           = malloc(NBASE * sizeof(double));
	PolyBivDFT* um_dft                = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBiv* um                       = malloc(poly_biv_bytes(params_glwe));
	double* um_univ_RnX               = calloc(NBASE, sizeof(double));

	// Define sk_ggsw = (1, 0, ... , 0)
	// TODO: why?
	sk_ggsw->values[0][0] = 1;

	// Computes the bivGGSW secret key out of the DFT domain
	transform_glwe_secret_key_not_dft_to_dft(module, sk_glwe_dft, sk_ggsw);

	// Draws uniformly both messages
	uniform_random_pol_znx(u_univ, NBASE, KAPPABASE);
	uniform_random_biv_poly(params_glwe_tilde, m, 1);

	//! Computation with function
	// Computes glwe_tilde, a bivGLWE(m) using the base-2Kappa_tilde decomposition
	glwe_secret_masking(module, glwe_tilde, sk_glwe_dft, m);

	// Computes ggsw, a bivGGSW(u) using the base-2Kappa
	glwegad_secret_encrypt(module, glwegad_ct, sk_glwe_dft, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition

	ggsw_external_product(module, ext_prod_computed, glwe_tilde, ggsw);
	normalize_glwe(module, ext_prod_computed, ext_prod_computed);

	// Computes the result phase = u*m + err , normalized with the base-2Kappa
	glwe_secret_demasking(module, phase_computed, sk_glwe_dft, ext_prod_computed);

	// The computed phase = u*m + err in Tn[X]
	biv_to_univ(params_glwe, um_computed_univ_RnX, phase_computed);

	//! Computation by hand
	// Computes DFT(m)
	pvda_vec_znx_dft(module, u_univ_dft, 1, u_univ, 1, NBASE);

	// Computes DFT(u*m)
	pvda_svp_apply_dft(module, um_dft, L_TILDEBASE, u_univ_dft, m, L_TILDEBASE, NBASE);

	// Computes u*m in ZN[X,Y]
	pvda_vec_znx_idft(module, um, L_TILDEBASE, um_dft, L_TILDEBASE);

	// Normalizes u*m with the base-2Kappa_tilde
	pvda_vec_znx_normalize_base2k(module, KAPPA_TILDEBASE, um, L_TILDEBASE, NBASE, um, L_TILDEBASE, NBASE);

	// Computes u*m in Tn[X]
	biv_to_univ(params_glwe_tilde, um_univ_RnX, um);

	//! Asserts um_computed_univ(X) = u * m_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(um_univ_RnX[p], um_computed_univ_RnX[p]);
		int cond    = diff < err_length;

		cr_assert(cond, "Equality failed with um_computed_univ_RnX[%ld] = %lf and  um_univ_RnX[%ld] = %lf", p, p,
		          um_univ_RnX[p] - floor(um_univ_RnX[p]), p, um_computed_univ_RnX[p], err_length);
	}

	// Clean up
	free(m);
	free(u_univ);
	free(u_univ_dft);
	free(phase_computed);
	free(um);
	free(um_univ_RnX);
	free(um_dft);
	free(um_computed_univ_RnX);

	delete_glwe(ext_prod_computed);
	delete_glwe(glwe_tilde);
	delete_glwegad(glwegad_ct);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_dft(sk_glwe_dft);
	delete_glwe_params(params_glwe);
	delete_glwe_params(params_glwe_tilde);
	delete_glwegad_params(params_glwegad);
	pvda_delete_module_info(module);
}
