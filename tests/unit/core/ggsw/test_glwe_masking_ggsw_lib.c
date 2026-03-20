#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "core/glwe/glwe_transform_key.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "rng.h"
#include "utils.h"

#define NBASE            1024
#define KBASE            8
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 4
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        -(LBASE / 2 + 1) * KAPPABASE

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 2
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)
#define SIGMA_TILDEBASE  -3

#define PROB_FACTOR      3

//! bivGGSW PART (begin)

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking_ggsw_lib, small_error)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// Since the message are drawn in Zn[X,Y], there is no decomposition error. Thus, the error should be smaller than 3*sigma 99.73% of the time
	double err_length = 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyBiv* m                           = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* err                         = malloc(poly_biv_bytes(params_glwe));
	PolyUnivRnX* m_univ_RnX              = calloc(NBASE, sizeof(double));
	PolyBiv* phase                       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	VecBiv* glwe_vec_computed            = malloc(glwe_coef_number(params_glwe) * sizeof(int64_t));
	PolyBiv* phase_computed              = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* phase_computed_univ_RnX = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// The input message uniformly drawn in Zn[X,Y]
	uniform_random_biv_poly(params_glwe, m, LBASE / 2);

	// The input error normaly drawn in Zn[X,Y]
	normal_random_biv_poly(params_glwe, err);

	//! Computations with functions
	// Computes the message in Tn[X] with the base-2Kappa
	biv_to_univ(params_glwe, m_univ_RnX, m);

	// The final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	GLWECiphertext glwe_ct = {params_glwe, glwe_vec_computed};
	glwe_secret_masking(module, &glwe_ct, sk_dft, phase);

	// The computed phase in Rn[X]
	glwe_secret_demasking(module, phase_computed, sk_dft, &glwe_ct);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Compare both phase in Rn[X]
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_computed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
	double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);

	//! Clean up
	free(phase_computed_univ_RnX);
	free(phase_computed);
	free(glwe_vec_computed);
	free(phase);
	free(m_univ_RnX);
	free(err);
	free(m);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking_ggsw_lib, uniform_RnX_message)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The message is drawn in Rn[X], there is a decomposition error of 2^{-kappa * l}. And the error should be smaller than 3*sigma 99.73% of the time.
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyUnivRnX* m_univ_RnX              = malloc(poly_univ_bytes(params_glwe));
	PolyBiv* err                         = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* m                           = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase                       = calloc(NBASE * LBASE, sizeof(int64_t));
	VecBiv* glwe_vec_computed            = malloc(glwe_coef_number(params_glwe) * sizeof(int64_t));
	PolyBiv* phase_computed              = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* phase_computed_univ_RnX = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Draws normaly in Rn[X] m_univ
	normal_random_vec(m_univ_RnX, NBASE, 0.0, 0.1);

	// Draws normaly the error
	normal_random_biv_poly(params_glwe, err);

	//! Computation with functions
	// Computes m_univ_RnXbivariate form
	univ_to_biv(params_glwe, m, m_univ_RnX);

	// Computes the final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the bivGLWE ciphertext
	GLWECiphertext glwe_ct = {params_glwe, glwe_vec_computed};
	glwe_secret_masking(module, &glwe_ct, sk_dft, phase);

	// Computes the computed phase in Rn[X]
	glwe_secret_demasking(module, phase_computed, sk_dft, &glwe_ct);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_computed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
	double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);

	//! Clean up
	free(phase_computed_univ_RnX);
	free(phase_computed);
	free(glwe_vec_computed);
	free(phase);
	free(m);
	free(err);
	free(m_univ_RnX);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

//! bivGGSW IN DFT SPACE PART (begin)

/**
 * @brief Test glwe_secret_masking_ggsw_lib_dft. In this test, the message is drawn in Zn[X,Y],
 * ie. there is no loss of precision for the message going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the message.
 */
Test(glwe_secret_masking_ggsw_lib_dft, small_error)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// Since the message are drawn in Zn[X,Y], there is no decomposition error. Thus, the error should be smaller than 3*sigma 99.73% of the time
	double err_length = 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyBiv* m                           = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* err                         = malloc(poly_biv_bytes(params_glwe));
	PolyUnivRnX* m_univ_RnX              = calloc(NBASE, sizeof(double));
	PolyBiv* phase                       = calloc(NBASE * LBASE, sizeof(int64_t));
	PolyBivDFT* phase_dft                = malloc(poly_biv_bytes(params_glwe));
	VecBivDFT* glwe_vec_computed_dft     = malloc(glwe_coef_number(params_glwe) * sizeof(int64_t));
	PolyBiv* phase_computed              = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* phase_computed_univ_RnX = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// The input message uniformly drawn in Zn[X,Y]
	uniform_random_biv_poly(params_glwe, m, LBASE / 2);

	// The input error normaly drawn in Zn[X,Y]
	normal_random_biv_poly(params_glwe, err);

	//! Computes with functions
	// Computes the message in Tn[X] with the base-2Kappa
	biv_to_univ(params_glwe, m_univ_RnX, m);

	// The final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the phase in the DFT domain
	pvda_vec_znx_dft(module, phase_dft, LBASE, phase, LBASE, NBASE);

	// Computes the bivGLWE ciphertext
	GLWECiphertextDFT glwe_dft_ct = {params_glwe, glwe_vec_computed_dft};
	glwe_secret_masking_dft(module, &glwe_dft_ct, sk_dft, phase_dft);

	// The computed phase in Rn[X]
	glwe_secret_demasking_dft(module, phase_computed, sk_dft, &glwe_dft_ct);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_computed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
	double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);

	//! Clean up
	free(phase_computed_univ_RnX);
	free(phase_computed);
	free(glwe_vec_computed_dft);
	free(phase_dft);
	free(phase);
	free(m_univ_RnX);
	free(err);
	free(m);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}

/**
 * @brief Test glwe_secret_masking_ggsw_lib. In this test, the message is drawn in Rn[X],
 * ie. there is a 2^(-l*kappa) loss of precision for the message, going from Rn[X] to Zn[X,Y].
 * Moreover, the error is small enough to not affect the bivariate message.
 */
Test(glwe_secret_masking_ggsw_lib_dft, uniform_RnX_message)
{
	// The variance of the error's distribution
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The message is drawn in Rn[X], there is a decomposition error of 2^{-kappa * l}. And the error should be smaller than 3*sigma 99.73% of the time.
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;

	//! Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);

	//! Variables
	GLWESecretKey* sk                    = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft             = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyUnivRnX* m_univ_RnX              = calloc(NBASE, sizeof(double));
	PolyBiv* err                         = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* m                           = malloc(poly_biv_bytes(params_glwe));
	PolyBiv* phase                       = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyBivDFT* phase_dft                = malloc(poly_biv_bytes(params_glwe));
	VecBivDFT* glwe_vec_computed_dft     = malloc(glwe_coef_number(params_glwe) * sizeof(int64_t));
	PolyBiv* phase_computed              = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* phase_computed_univ_RnX = calloc(NBASE, sizeof(double));

	//! Draws each input variable
	// Draws uniformly in (Cm[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Draws normaly in Rn[X] m_univ_RnX
	normal_random_vec(m_univ_RnX, NBASE, 0.0, 0.1);

	// Draws normaly the error
	normal_random_biv_poly(params_glwe, err);

	//! Computation with functions
	// Computes m_univ_RnX bivariate form
	univ_to_biv(params_glwe, m, m_univ_RnX);

	// Computes the final phase = m + err
	add_biv_poly(params_glwe, phase, NBASE, m, NBASE, err, NBASE);

	// Computes the phase in the DFT domain
	pvda_vec_znx_dft(module, phase_dft, LBASE, phase, LBASE, NBASE);

	// Computes the bivGLWE ciphertext
	GLWECiphertextDFT glwe_dft_ct = {params_glwe, glwe_vec_computed_dft};
	glwe_secret_masking_dft(module, &glwe_dft_ct, sk_dft, phase_dft);

	// Computes the computed phase in Rn[X]
	glwe_secret_demasking_dft(module, phase_computed, sk_dft, &glwe_dft_ct);

	// The computed phase in Rn[X]
	biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

	// A variable counting the number of times the error is greater than 3*sigma
	int big_error_count = 0;

	// Using the triangle inequality, for each p, the difference should be smaller than |err_p| + |msg_p -
	// msgComputed_p| Ie, 3*sigma + 2^(-l*kappa)
	for (uint64_t p = 0; p < NBASE; p++)
	{
		double diff = torus_distance(m_univ_RnX[p], phase_computed_univ_RnX[p]);

		int cond = diff < err_length;

		if (!cond) big_error_count++;
	}

	int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
	double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

	/// Asserts big_error_count <= 0.0027*N
	cr_assert(big_error_count <= max_fails,
	          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
	          "chance, that happens.",
	          max_fails, big_error_count, proba);

	//! Clean up
	free(phase_computed_univ_RnX);
	free(phase_computed);
	free(glwe_vec_computed_dft);
	free(phase_dft);
	free(phase);
	free(m);
	free(err);
	free(m_univ_RnX);
	delete_glwe_secret_key(sk);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	pvda_delete_module_info(module);
}
