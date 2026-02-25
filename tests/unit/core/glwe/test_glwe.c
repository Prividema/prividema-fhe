#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "common/spqlios_alias.h"
#include "core/glwe/glwe.h"
#include "rng.h"

#define NBASE      1024
#define KBASE      1
#define KAPPABASE  4
#define NLIMBSBASE (KBASE + 1) * 8
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABASE  -(LBASE / 2 + 1) * KAPPABASE

//! GGWS PART (begin)

/**
 * @brief Test glwe_secret_masking. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking, small_error)
{
	// The variance of the error's distribution
	double sigma             = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	GLWECtParams* params     = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	MODULE* module           = new_module_info_p(NBASE);

	GLWECiphertext* ct       = new_glwe(params);
	GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(module, KBASE, 3);

	// The input message uniformly drawn in Zn[X,Y]
	PolyBiv* msg     = new_uniform_random_biv_poly(module, params, LBASE / 2);

	double* msg_univ = calloc(NBASE, sizeof(double));
	biv_to_univ(params, msg_univ, msg);

	// The input error normaly drawn in Zn[X,Y]
	PolyBiv* err = new_normal_random_biv_poly(module, params);

	// The final phase = msg + err
	PolyBiv* phase = calloc(NBASE * LBASE, sizeof(int64_t));
	add_biv_poly(params, phase, NBASE, msg, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking(module, ct, sk_dft, phase);

	// The computed phase in Rn[X]
	PolyBiv* phase_computed = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	glwe_secret_demasking(module, phase_computed, sk_dft, ct);

	// The computed phase in Rn[X]
	double* phase_computed_univ = calloc(NBASE, sizeof(double));
	biv_to_univ(params, phase_computed_univ, phase_computed);

	// Compare both phase in Rn[X]
	for (int64_t p = 0; p < NBASE; p++) {
		double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
		double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma;

		int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
		cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p,
		          msg_univ[p], p, phase_computed_univ[p]);
	}

	free(msg);
	free(msg_univ);
	free(err);
	free(phase);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe(ct);
	delete_module_info_p(module);
	delete_glwe_ct_params(params);
	delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Test glwe_encrypt_priv. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking, uniform_RnX_message)
{
	// The variance of the error's distribution
	double sigma             = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	GLWECtParams* params     = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	MODULE* module           = new_module_info_p(NBASE);

	GLWECiphertext* ct       = new_glwe(params);
	GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(module, KBASE, 3);

	// The input message uniformly drawn in Rn[X]
	double* msg_univ = malloc(poly_univ_bytes(params));
	new_normal_random_vec(NBASE, msg_univ, 1, NBASE, 0.0, 0.1);

	// Computes the bivariate form
	PolyBiv* msg = malloc(poly_biv_bytes(params));
	univ_to_biv(params, msg, msg_univ);

	// The input error normaly drawn in Zn[X,Y]
	PolyBiv* err = new_normal_random_biv_poly(module, params);

	// The final phase = msg + err
	PolyBiv* phase = calloc(NBASE * LBASE, sizeof(int64_t));
	add_biv_poly(params, phase, NBASE, msg, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking(module, ct, sk_dft, phase);

	// The computed phase in Rn[X]
	PolyBiv* phase_computed = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	glwe_secret_demasking(module, phase_computed, sk_dft, ct);

	// The computed phase in Rn[X]
	double* phase_computed_univ = calloc(NBASE, sizeof(double));
	biv_to_univ(params, phase_computed_univ, phase_computed);

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, then |err_p| + 2^(-l*kappa) Assures the error, of length (-lN/2), affects the message, of
	// degree l/2 in Y
	for (int64_t p = 0; p < NBASE; p++) {
		double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
		double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma;

		int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
		cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p,
		          msg_univ[p], p, phase_computed_univ[p]);
	}

	free(msg);
	free(msg_univ);
	free(err);
	free(phase);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe(ct);
	delete_module_info_p(module);
	delete_glwe_ct_params(params);
	delete_glwe_secret_key_dft(sk_dft);
}

//! GLWE IN DFT PART (begin)
/**
 * @brief Test glwe_secret_masking. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking_dft, small_error)
{
	// The variance of the error's distribution
	double sigma              = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	GLWECtParams* params      = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	MODULE* module            = new_module_info_p(NBASE);

	GLWECiphertextDFT* ct_dft = new_glwe_dft(params);
	GLWESecretKeyDFT* sk_dft  = new_uniform_glwe_secret_key_dft(module, KBASE, 3);

	// The input message uniformly drawn in Zn[X,Y]
	PolyBiv* msg     = new_uniform_random_biv_poly(module, params, LBASE / 2);

	double* msg_univ = calloc(NBASE, sizeof(double));
	biv_to_univ(params, msg_univ, msg);

	// The input error normaly drawn in Zn[X,Y]
	PolyBiv* err = new_normal_random_biv_poly(module, params);

	// The final phase = msg + err
	PolyBiv* phase = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	add_biv_poly(params, phase, NBASE, msg, NBASE, err, NBASE);

	// The phase in DFT space
	PolyBivDFT* phase_dft = calloc(poly_biv_coef_number(params), sizeof(double));
	vec_znx_dft_p(module, phase_dft, LBASE, phase, LBASE, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking_dft(module, ct_dft, sk_dft, phase_dft);

	// The computed phase in Rn[X]
	PolyBiv* phase_computed = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	glwe_secret_demasking_dft(module, phase_computed, sk_dft, ct_dft);

	// The computed phase in Rn[X]
	double* phase_computed_univ = calloc(NBASE, sizeof(double));
	biv_to_univ(params, phase_computed_univ, phase_computed);

	// Assures the error, of length (-lN/2), affects the message, of degree l/2 in Y
	for (int64_t p = 0; p < NBASE; p++) {
		double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
		double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma;

		int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
		cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p,
		          msg_univ[p], p, phase_computed_univ[p]);
	}

	free(msg);
	free(msg_univ);
	free(err);
	free(phase);
	free(phase_dft);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe_dft(ct_dft);
	delete_module_info_p(module);
	delete_glwe_ct_params(params);
	delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Test glwe_secret_masking_dft. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking_dft, uniform_RnX_message)
{
	// The variance of the error's distribution
	double sigma              = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	GLWECtParams* params      = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	MODULE* module            = new_module_info_p(NBASE);

	GLWECiphertextDFT* ct_dft = new_glwe_dft(params);
	GLWESecretKeyDFT* sk_dft  = new_uniform_glwe_secret_key_dft(module, KBASE, 3);

	// The input message uniformly drawn in Rn[X]
	double* msg_univ = malloc(poly_univ_bytes(params));
	new_normal_random_vec(NBASE, msg_univ, 1, NBASE, 0.0, 0.1);

	// Computes the bivariate form
	PolyBiv* msg = malloc(poly_biv_bytes(params));
	univ_to_biv(params, msg, msg_univ);

	// The input error normaly drawn in Zn[X,Y]
	PolyBiv* err = new_normal_random_biv_poly(module, params);

	// The final phase = msg + err
	PolyBiv* phase = calloc(NBASE * LBASE, sizeof(int64_t));
	add_biv_poly(params, phase, NBASE, msg, NBASE, err, NBASE);

	// The phase in DFT space
	PolyBivDFT* phase_dft = malloc(poly_biv_bytes(params));
	vec_znx_dft_p(module, phase_dft, LBASE, phase, LBASE, NBASE);

	// Computes the bivGLWE ciphertext
	glwe_secret_masking_dft(module, ct_dft, sk_dft, phase_dft);

	// The computed phase in Rn[X]
	PolyBiv* phase_computed = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	glwe_secret_demasking_dft(module, phase_computed, sk_dft, ct_dft);

	// The computed phase in Rn[X]
	double* phase_computed_univ = calloc(NBASE, sizeof(double));
	biv_to_univ(params, phase_computed_univ, phase_computed);

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, then |err_p| + 2^(-l*kappa)
	for (int64_t p = 0; p < NBASE; p++) {
		double diff_1 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p];
		double diff_2 = msg_univ[p] - round(msg_univ[p]) - phase_computed_univ[p] + floor(phase_computed_univ[p]) +
		                ceil(phase_computed_univ[p]);
		double err_length = 3 * sigma + ldexp(1.0, -LBASE * KAPPABASE);

		int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);
		cr_assert(cond, "Equality failed at p = %ld with msg_univ[%ld] = %lf and phase_computed_univ[%ld] = %lf", p, p,
		          msg_univ[p], p, phase_computed_univ[p]);
	}

	free(msg);
	free(msg_univ);
	free(err);
	free(phase);
	free(phase_dft);
	free(phase_computed);
	free(phase_computed_univ);
	delete_glwe_dft(ct_dft);
	delete_module_info_p(module);
	delete_glwe_ct_params(params);
	delete_glwe_secret_key_dft(sk_dft);
}
