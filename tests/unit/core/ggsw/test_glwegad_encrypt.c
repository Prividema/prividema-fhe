#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw.h"
#include "core/ggsw/glwegad.h"
#include "core/glwe/glwe.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_transform_key.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

#define NBASE            1024
#define KBASE            1
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 4
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        -(LBASE / 2 + 1) * KAPPABASE

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 4
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)
#define SIGMA_TILDEBASE  -3

#define PROB_FACTOR      3

/**
 * @brief Tests ggsw_secret_encrpyt
 *
 */
Test(glwegad_secret_encrypt, works)
{
	// bivGLWE and bivGGSW parameters. This set of bivGLWE parameters is for bivGGSW ciphertext
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99%
	// chance
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;
	cr_log_info("error length = %e", err_length);

	// Parameters
	MODULE* module                   = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe          = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GLWEGadgetParams* params_glwegad = new_glwegad_params(params_glwe, KAPPA_TILDEBASE, L_TILDEBASE);

	// Variables
	GLWESecretKey* sk             = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft      = alloc_glwe_secret_key_dft(NBASE, KBASE);
	GLWEGadgetCiphertext* glwegad = new_glwegad(params_glwegad);
	PolyUniv* m_univ              = new_univ(params_glwe);
	PolyUnivDFT* m_univ_dft       = new_univ_dft(module);
	// Variables to compute the phase of each ggsw's row
	PolyBiv* phase_computed              = new_biv_poly(params_glwe);
	PolyUnivRnX* phase_computed_univ_RnX = new_univ_rnx(params_glwe);
	PolyUnivRnX* phase_univ_RnX          = new_univ_rnx(params_glwe);
	PolyUnivDFT* m_skj_univ_dft          = new_univ_dft(module);
	PolyUniv* m_skj_univ                 = new_univ(params_glwe);

	// Draws uniformly in (Zn[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Draws uniformly in Zn[X] the message
	uniform_random_pol_znx(m_univ, NBASE, KAPPABASE);

	// Computes a bivGGSW(m)
	glwegadget_secret_encrypt(module, glwegad, sk_dft, m_univ);

	for (uint64_t i = 1; i <= L_TILDEBASE; i++)
	{
		// Fills each changed variable with 0s'
		memset(phase_computed, 0, poly_biv_bytes(params_glwe));
		memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
		memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

		// The bivGGSW row correponding to bivGLWE(m / 2^{kappa_tilde * i})

		// Point to bivGLWE(m / 2^{kappa_tilde * i})

		VecBiv* glwe_vec_ptr = glwegad->mat + (i - 1) * glwe_coef_number(params_glwe);

		// Computes the phase = m / 2^{kappa_tilde * i} + err
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec_ptr};
		glwe_secret_demasking(module, phase_computed, sk_dft, &glwe_ct);

		// Computes the phase in
		biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

		//! Computes by hand the phase = m / 2^{kappa_tilde * i}
		for (uint64_t p = 0; p < NBASE; p++)
		{
			phase_univ_RnX[p] = ldexp((double)m_univ[p], -(params_glwegad->kappa_tilde * i));
		}

		// A variable counting the number of times the error is greater than 3*sigma
		int big_error_count = 0;

		// Assures that the difference between the phase = m / 2^{kappa_tilde * i} and the computed phase,
		// are only different by an error of approximation and a gaussian error
		for (uint64_t p = 0; p < NBASE; p++)
		{
			double diff = torus_distance(phase_univ_RnX[p], phase_computed_univ_RnX[p]);

			int cond = diff < err_length;

			if (!cond) big_error_count++;
		}

		/// Prob that the number of error grater than 3sigma is greater or equal than 0.0027*N
		int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
		double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

		/// Asserts big_error_count <= 0.0027*N
		cr_assert(big_error_count <= max_fails,
		          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
		          "chance, that happens.",
		          max_fails, big_error_count, proba);
	}

	// Clean up
	delete_univ_rnx(phase_computed_univ_RnX);
	delete_univ_dft(m_skj_univ_dft);
	delete_univ(m_skj_univ);
	free(phase_computed);
	delete_univ_rnx(phase_univ_RnX);
	delete_univ(m_univ);
	delete_univ_dft(m_univ_dft);
	delete_glwegad(glwegad);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	delete_glwegad_params(params_glwegad);
	pvda_delete_module_info(module);
}
