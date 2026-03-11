#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ciphertext.h"
#include "rng.h"

#define NBASE            1024
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

//! GGSW PART (begin)

/** The test is done without error, it is a proof of concept*/
Test(ggsw_external_product, without_error)
{
	//! Variance of the error's normal distributions
	double sigma       = 0;
	double sigma_tilde = 0;

	//! Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GLWECtParams* params_glwe_tilde =
	    new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
	GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
	MODULE* module            = new_module_info(NBASE, FFT64);

	//! Variables
	GGSWSecretKey* sk_ggsw            = new_ggsw_secret_key(NBASE, KBASE);
	GGSWSecretKeyDFT* sk_ggsw_dft     = new_ggsw_secret_key_dft(NBASE, KBASE);
	GLWESecretKeyDFT* sk_glwe_dft     = new_glwe_secret_key_dft(NBASE, KBASE);
	GGSWCiphertext* ggsw              = new_ggsw(params_ggsw);
	GLWECiphertext* glwe_tilde        = new_glwe(params_glwe_tilde);
	GLWECiphertext* ext_prod_computed = new_glwe(params_glwe);
	PolyUniv* u_univ                  = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* m                        = malloc(poly_biv_bytes(params_glwe_tilde));

	//! Variables to compute the result phase of the external product
	PolyBiv* phase_computed       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* um_univ_computed = calloc(NBASE, sizeof(double));
	PolyUnivDFT* u_univ_dft       = malloc(NBASE * sizeof(double));
	PolyBivDFT* um_dft            = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBiv* um                   = malloc(poly_biv_bytes(params_glwe));
	double* um_univ               = calloc(NBASE, sizeof(double));

	// Define sk_ggsw = (1, 0, ... , 0)
	sk_ggsw->values[0][0] = 1;
	sk_ggsw->values[0][1] = 0;

	// Computes the GGSW secret key out of the DFT domain
	transform_ggsw_secret_key_not_dft_to_dft(module, sk_ggsw_dft, sk_ggsw);
	transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_glwe_dft, sk_ggsw_dft);

	// Draws uniformly both messages
	uniform_random_pol_znx(u_univ, NBASE, KAPPABASE);
	uniform_random_biv_poly(module, params_glwe_tilde, m, 1);

	//! Computation with function
	// Computes glwe_tilde, a bivGLWE(m) using the base-2Kappa_tilde decomposition
	glwe_secret_masking(module, glwe_tilde, sk_glwe_dft, m);

	// Computes ggsw, a bivGGSW(u) using the base-2Kappa
	ggsw_secret_encrypt(module, params_ggsw, ggsw, sk_ggsw_dft, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
	ggsw_external_product(module, ext_prod_computed, glwe_tilde, ggsw);
	normalize_glwe(module, ext_prod_computed, ext_prod_computed);

	// Computes the result phase = u*m + err , normalized with the base-2Kappa
	glwe_secret_demasking(module, phase_computed, sk_glwe_dft, ext_prod_computed);

	// The computed phase = u*m + err in Tn[X]
	biv_to_univ(params_glwe, um_univ_computed, phase_computed);

	//! Computation by hand
	// Computes DFT(m)
	vec_znx_dft_p(module, u_univ_dft, 1, u_univ, 1, NBASE);

	// Computes DFT(u*m)
	svp_apply_dft_p(module, um_dft, L_TILDEBASE, u_univ_dft, m, L_TILDEBASE, NBASE);

	// Computes u*m in ZN[X,Y]
	vec_znx_idft_p(module, um, L_TILDEBASE, um_dft, L_TILDEBASE);

	// Normalizes u*m with the base-2Kappa_tilde
	vec_znx_normalize_base2k_p(module, KAPPA_TILDEBASE, um, L_TILDEBASE, NBASE, um, L_TILDEBASE, NBASE);

	// Computes u*m in Tn[X]
	biv_to_univ(params_glwe_tilde, um_univ, um);

	//! Asserts um_univ_computed(X) = u * m_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff_1 = um_univ[p] - floor(um_univ[p]) - um_univ_computed[p];
		double diff_2 = um_univ[p] - floor(um_univ[p]) - um_univ_computed[p] + floor(um_univ_computed[p]) +
		                ceil(um_univ_computed[p]);
		double err_length = ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE);

		int cond = (diff_1 <= err_length && diff_1 >= -err_length) || (diff_2 <= err_length && diff_2 >= -err_length);

		cr_assert(cond, "Equality failed with um_univ_computed[%ld] = %lf and  um_univ[%ld] = %lf", p, p,
		          um_univ[p] - floor(um_univ[p]), p, um_univ_computed[p] - floor(um_univ_computed[p]), err_length);
	}

	// Clean up
	free(m);
	free(u_univ);
	free(u_univ_dft);
	free(phase_computed);
	free(um);
	free(um_univ);
	free(um_dft);
	free(um_univ_computed);

	delete_glwe(ext_prod_computed);
	delete_glwe(glwe_tilde);
	delete_ggsw(ggsw);

	delete_ggsw_secret_key(sk_ggsw);
	delete_glwe_secret_key_dft(sk_glwe_dft);
	delete_ggsw_secret_key_dft(sk_ggsw_dft);

	delete_glwe_ct_params(params_glwe);
	delete_glwe_ct_params(params_glwe_tilde);
	delete_ggsw_ct_params(params_ggsw);
	delete_module_info_p(module);
}

//! GGSW PART in the DFT domain (begin)

Test(ggsw_external_product_dft, without_error)
{
	//! Variance of the error's normal distributions
	double sigma       = 0;
	double sigma_tilde = 0;

	//! Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GLWECtParams* params_glwe_tilde =
	    new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
	GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
	MODULE* module            = new_module_info(NBASE, FFT64);

	//! Variables
	GGSWSecretKey* sk_ggsw                = new_ggsw_secret_key(NBASE, KBASE);
	GGSWSecretKeyDFT* sk_ggsw_dft         = new_ggsw_secret_key_dft(NBASE, KBASE);
	GLWESecretKeyDFT* sk_glwe_dft         = new_glwe_secret_key_dft(NBASE, KBASE);
	GGSWCiphertext* ggsw_dft              = new_ggsw(params_ggsw);
	GLWECiphertext* glwe_tilde_dft        = new_glwe(params_glwe_tilde);
	GLWECiphertext* ext_prod_computed_dft = new_glwe(params_glwe);
	PolyUniv* u_univ                      = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* m                            = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBivDFT* m_dft                     = malloc(poly_biv_bytes(params_glwe));

	//! Variables to compute the result phase of the external product
	PolyBiv* phase_computed       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* um_univ_computed = calloc(NBASE, sizeof(double));
	PolyUnivDFT* u_univ_dft       = malloc(NBASE * sizeof(double));
	PolyBivDFT* um_dft            = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBiv* um                   = malloc(poly_biv_bytes(params_glwe));
	double* um_univ               = calloc(NBASE, sizeof(double));

	// Define sk_ggsw = (1, 0, ... , 0)
	sk_ggsw->values[0][0] = 1;
	sk_ggsw->values[0][1] = 0;

	// Computes the GGSW secret key out of the DFT domain
	transform_ggsw_secret_key_not_dft_to_dft(module, sk_ggsw_dft, sk_ggsw);
	transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(sk_glwe_dft, sk_ggsw_dft);

	// Draws uniformly both messages
	uniform_random_pol_znx(u_univ, NBASE, KAPPABASE);
	uniform_random_biv_poly(module, params_glwe_tilde, m, 1);

	//! Computation with function
	// Computes glwe_tilde, a bivGLWE(m) using the base-2Kappa_tilde decomposition
	glwe_secret_masking_dft(module, glwe_tilde_dft, sk_glwe_dft, m_dft);

	// Computes ggsw, a bivGGSW(u) using the base-2Kappa
	ggsw_secret_encrypt_dft(module, params_ggsw, ggsw_dft, sk_ggsw_dft, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
	ggsw_external_product_dft(module, ext_prod_computed_dft, glwe_tilde_dft, ggsw_dft);
	normalize_glwe_dft(module, ext_prod_computed_dft, ext_prod_computed_dft);

	// Computes the result phase = u*m + err , normalized with the base-2Kappa
	glwe_secret_demasking_dft(module, phase_computed, sk_glwe_dft, ext_prod_computed_dft);

	// The computed phase = u*m + err in Tn[X]
	biv_to_univ(params_glwe, um_univ_computed, phase_computed);

	//! Computation by hand
	// Computes DFT(m)
	vec_znx_dft_p(module, u_univ_dft, 1, u_univ, 1, NBASE);

	// Computes DFT(u*m)
	svp_apply_dft_p(module, um_dft, L_TILDEBASE, u_univ_dft, m, L_TILDEBASE, NBASE);

	// Computes u*m in ZN[X,Y]
	vec_znx_idft_p(module, um, L_TILDEBASE, um_dft, L_TILDEBASE);

	// Normalizes u*m with the base-2Kappa_tilde
	vec_znx_normalize_base2k_p(module, KAPPA_TILDEBASE, um, L_TILDEBASE, NBASE, um, L_TILDEBASE, NBASE);

	// Computes u*m in Tn[X]
	biv_to_univ(params_glwe_tilde, um_univ, um);

	//! Asserts um_univ_computed(X) = u * m_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff_1 = um_univ[p] - floor(um_univ[p]) - um_univ_computed[p];
		double diff_2 = um_univ[p] - floor(um_univ[p]) - um_univ_computed[p] + floor(um_univ_computed[p]) +
		                ceil(um_univ_computed[p]);
		double err_length = ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE);

		int cond = (diff_1 <= err_length && diff_1 >= -err_length) || (diff_2 <= err_length && diff_2 >= -err_length);

		cr_assert(cond, "Equality failed with um_univ_computed[%ld] = %lf and  um_univ[%ld] = %lf", p, p,
		          um_univ[p] - floor(um_univ[p]), p, um_univ_computed[p] - floor(um_univ_computed[p]), err_length);
	}

	// Clean up
	free(m);
	free(u_univ);
	free(u_univ_dft);
	free(phase_computed);
	free(um);
	free(um_univ);
	free(um_dft);
	free(um_univ_computed);

	delete_glwe_dft(ext_prod_computed_dft);
	delete_glwe_dft(glwe_tilde_dft);
	delete_ggsw_dft(ggsw_dft);

	delete_ggsw_secret_key(sk_ggsw);
	delete_glwe_secret_key_dft(sk_glwe_dft);
	delete_ggsw_secret_key_dft(sk_ggsw_dft);

	delete_glwe_ct_params(params_glwe);
	delete_glwe_ct_params(params_glwe_tilde);
	delete_ggsw_ct_params(params_ggsw);
	delete_module_info_p(module);
}