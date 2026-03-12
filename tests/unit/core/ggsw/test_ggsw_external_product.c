#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ciphertext.h"
#include "core/glwe/glwe_transform_key.h"
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

//! bivGGSW PART (begin)

void printf_glwe(GLWECiphertext* glwe){
	for(int64_t j = 0 ; j < KBASE+1 ; j++ )
		printf_poly_biv(glwe->vec+ j*NBASE, (KBASE + 1) * NBASE, NBASE, LBASE);
}

void printf_glwe_dft(MODULE* module, GLWECiphertextDFT* glwe_dft)
{
	VecBiv* glwe_vec = malloc(glwe_bytes(glwe_dft->params));
	vec_znx_idft_p(module, glwe_vec, glwe_size(glwe_dft->params), glwe_dft->vec, glwe_size(glwe_dft->params));

	for(int64_t j = 0 ; j < KBASE+1 ; j++ )
		printf_poly_biv(glwe_vec + j*NBASE, (KBASE + 1) * NBASE, NBASE, LBASE);

	free(glwe_vec);
}

/** The test is done without error, it is a proof of concept*/
Test(ggsw_external_product, without_error)
{
	//! Variance of the error's normal distributions
	double sigma       = 0;
	double sigma_tilde = 0;
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;

	//! Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GLWECtParams* params_glwe_tilde =
	    new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
	GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
	MODULE* module            = new_module_info(NBASE, FFT64);

	//! Variables
	GLWESecretKey* sk_ggsw            = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_glwe_dft     = alloc_glwe_secret_key_dft(NBASE, KBASE);
	GGSWCiphertext* ggsw              = new_ggsw(params_ggsw);
	GLWECiphertext* glwe_tilde        = new_glwe(params_glwe_tilde);
	GLWECiphertext* ext_prod_computed = new_glwe(params_glwe);
	PolyUniv* u_univ                  = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* m                        = malloc(poly_biv_bytes(params_glwe_tilde));

	//! Variables to compute the result phase of the external product
	PolyBiv* phase_computed       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* um_computed_univ_RnX = calloc(NBASE, sizeof(double));
	PolyUnivDFT* u_univ_dft       = malloc(NBASE * sizeof(double));
	PolyBivDFT* um_dft            = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBiv* um                   = malloc(poly_biv_bytes(params_glwe));
	double* um_univ_RnX               = calloc(NBASE, sizeof(double));

	// Define sk_ggsw = (1, 0, ... , 0)
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
	ggsw_secret_encrypt(module, params_ggsw, ggsw, sk_glwe_dft, u_univ);

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
	vec_znx_dft_p(module, u_univ_dft, 1, u_univ, 1, NBASE);

	// Computes DFT(u*m)
	svp_apply_dft_p(module, um_dft, L_TILDEBASE, u_univ_dft, m, L_TILDEBASE, NBASE);

	// Computes u*m in ZN[X,Y]
	vec_znx_idft_p(module, um, L_TILDEBASE, um_dft, L_TILDEBASE);

	// Normalizes u*m with the base-2Kappa_tilde
	vec_znx_normalize_base2k_p(module, KAPPA_TILDEBASE, um, L_TILDEBASE, NBASE, um, L_TILDEBASE, NBASE);

	// Computes u*m in Tn[X]
	biv_to_univ(params_glwe_tilde, um_univ_RnX, um);

	//! Asserts um_computed_univ(X) = u * m_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{

    double diff = torus_distance(um_univ_RnX[p], um_computed_univ_RnX[p]);
    int cond = diff < err_length;

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
	delete_ggsw(ggsw);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_dft(sk_glwe_dft);
	delete_glwe_ct_params(params_glwe);
	delete_glwe_ct_params(params_glwe_tilde);
	delete_ggsw_ct_params(params_ggsw);
	delete_module_info_p(module);
}

//! bivGGSW PART in the DFT domain (begin)

Test(ggsw_external_product_dft, without_error)
{
	//! Variance of the error's normal distributions
	double sigma       = 0;
	double sigma_tilde = 0;
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;

	//! Parameters
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GLWECtParams* params_glwe_tilde =
	    new_glwe_ct_params(NBASE, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE, SIGMA_TILDEBASE);
	GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
	MODULE* module            = new_module_info(NBASE, FFT64);

	//! Variables
	GLWESecretKey* sk_ggsw                = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_glwe_dft         = alloc_glwe_secret_key_dft(NBASE, KBASE);
	GGSWCiphertextDFT* ggsw_dft              = new_ggsw_dft(params_ggsw);
	GLWECiphertextDFT* glwe_tilde_dft        = new_glwe_dft(params_glwe_tilde);
	GLWECiphertextDFT* ext_prod_computed_dft = new_glwe_dft(params_glwe);
	PolyUniv* u_univ                     = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* m                            = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBivDFT* m_dft                     = malloc(poly_biv_bytes(params_glwe));

	//! Variables to compute the result phase of the external product
	PolyBiv* phase_computed       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* um_computed_univ_RnX = calloc(NBASE, sizeof(double));
	PolyUnivDFT* u_univ_dft       = malloc(NBASE * sizeof(double));
	PolyBivDFT* um_dft            = malloc(poly_biv_bytes(params_glwe_tilde));
	PolyBiv* um                   = malloc(poly_biv_bytes(params_glwe));
	PolyUnivRnX* um_univ_RnX               = calloc(NBASE, sizeof(double));

	// Define sk_ggsw = (1, 0, ... , 0)
	sk_ggsw->values[0][0] = 2;

	// Computes the bivGGSW secret key out of the DFT domain
	transform_glwe_secret_key_not_dft_to_dft(module, sk_glwe_dft, sk_ggsw);

	// Draws uniformly both messages
	uniform_random_pol_znx(u_univ, NBASE, KAPPABASE);
	uniform_random_biv_poly(params_glwe_tilde, m, 1);

	// Computes m in the DFT space
	vec_znx_dft_p(module, m_dft, 1, m, 1, NBASE);

	//! Computation with function
	// Computes glwe_tilde, a bivGLWE(m) using the base-2Kappa_tilde decomposition
	glwe_secret_masking_dft(module, glwe_tilde_dft, sk_glwe_dft, m_dft);

	// Computes ggsw, a bivGGSW(u) using the base-2Kappa
	ggsw_secret_encrypt_dft(module, params_ggsw, ggsw_dft, sk_glwe_dft, u_univ);

	// Computes the external product of glwe_tilde and ggsw
	// It should result in a bivGLWE(u*m) using the base-2Kappa decomposition
	ggsw_external_product_dft(module, ext_prod_computed_dft, glwe_tilde_dft, ggsw_dft);
	normalize_glwe_dft(module, ext_prod_computed_dft, ext_prod_computed_dft);

	// Computes the result phase = u*m + err , normalized with the base-2Kappa
	glwe_secret_demasking_dft(module, phase_computed, sk_glwe_dft, ext_prod_computed_dft);

	// The computed phase = u*m + err in Tn[X]
	biv_to_univ(params_glwe, um_computed_univ_RnX, phase_computed);

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
	biv_to_univ(params_glwe_tilde, um_univ_RnX, um);

	//! Asserts um_computed_univ(X) = u * m_univ
	for (uint64_t p = 0; p < NBASE; p++)
	{

    double diff = torus_distance(um_univ_RnX[p], um_computed_univ_RnX[p]);
    int cond = diff < err_length;

		cr_assert(cond, "Equality failed with um_computed_univ_RnX[%ld] = %lf and  um_univ_RnX[%ld] = %lf", p, p,
		          um_computed_univ_RnX[p], p, um_univ_RnX[p], err_length);
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

	delete_glwe_dft(ext_prod_computed_dft);
	delete_glwe_dft(glwe_tilde_dft);
	delete_ggsw_dft(ggsw_dft);

	delete_glwe_secret_key(sk_ggsw);
	delete_glwe_secret_key_dft(sk_glwe_dft);

	delete_glwe_ct_params(params_glwe);
	delete_glwe_ct_params(params_glwe_tilde);
	delete_ggsw_ct_params(params_ggsw);
	delete_module_info_p(module);
}
