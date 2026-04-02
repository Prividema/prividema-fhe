#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "bivariate_polynomial.h"
#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "ggsw_ciphertext.h"
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
Test(ggsw_secret_encrypt, works)
{
	// bivGLWE and bivGGSW parameters. This set of bivGLWE parameters is for bivGGSW ciphertext
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99%
	// chance
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;
	cr_log_info("error length = %e", err_length);

	// Parameters
	MODULE* module          = pvda_new_module_info(NBASE);
	GLWEParams* params_glwe = new_glwe_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GGSWParams* params_ggsw = new_ggsw_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GLWESecretKey* sk        = alloc_glwe_secret_key(params_glwe);
	GLWESecretKeyDFT* sk_dft = alloc_glwe_secret_key_dft(params_glwe);
	GGSWCiphertext* ggsw     = new_ggsw(params_ggsw);
	PolyUniv* m_univ         = new_univ(params_glwe);
	PolyUnivDFT* m_univ_dft  = new_univ_dft(module);
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

	// Computes the message in the DFT domain
	univ_coefs_to_dft(module, m_univ_dft, m_univ);

	// Computes a bivGGSW(m)
	ggsw_secret_encrypt(module, ggsw, sk_dft, m_univ);

	// Asserts the j-th row in the i-th PartialGGSW(m) in the ggsw is :
	// - a bivGLWE(-m * sk_j / 2^{kappa_tilde * i})), for j < k
	// - a bivGLWE(m / 2^{kappa_tilde * i}))        , for j = k
	for (uint64_t i = 1; i <= L_TILDEBASE; i++)
	{
		// For j from 0 to k-1 (ie, the -m*sk values)
		for (uint64_t j = 0; j < K_TILDEBASE; j++)
		{
			// Fills each changed variable with 0s'
			memset(phase_computed, 0, poly_biv_bytes(params_glwe));
			memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ_dft, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ, 0, poly_univ_bytes(params_glwe));
			memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

			// Computes the phase = -m * sk_j / 2^{kappa_tilde * i}) + err
			GLWECiphertext glwe_ct = {params_glwe, ggsw_retrieve_bivglwe(ggsw, j, i)};
			glwe_secret_decrypt(module, phase_computed, sk_dft, &glwe_ct);
			biv_to_univ_rnx(params_glwe, phase_computed_univ_RnX, phase_computed);

			// Computes DFT(m * sk_j)
			mult_vec_znx_dft(module, m_skj_univ_dft, 1, glwe_sk_extract_poly_dft(sk_dft, j), 1, m_univ_dft, 1);

			// Computes DFT(-m * sk_j)
			// TODO: znx negate
			for (uint64_t p = 0; p < NBASE; p++)
			{
				m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
			}
			pvda_vec_znx_idft(module, m_skj_univ, 1, m_skj_univ_dft, 1);

			// Computes -m * sk_j / 2^{kappa_tilde * i}
			for (uint64_t p = 0; p < NBASE; p++)
			{
				phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -(params_ggsw->kappa_tilde * i));
			}

			// A variable counting the number of times the error is greater than 3*sigma
			int big_error_count = 0;

			// Assures that the difference between the phase = m / 2^{kappa_tilde * i} and the computed phase,
			// are only different by an error of approximation and a gaussian error
			for (uint64_t p = 0; p < NBASE; p++)
			{
				double distance = torus_distance(phase_univ_RnX[p], phase_computed_univ_RnX[p]);
				int cond        = distance < err_length;

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

		// Final row (k/k+1, ie, the b in (a0, a1, ..., ak-1, b))

		// Fills each changed variable with 0s'
		memset(phase_computed, 0, poly_biv_bytes(params_glwe));
		memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
		memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

		// The bivGGSW row correponding to bivGLWE(m / 2^{kappa_tilde * i})
		uint64_t row_i_ktilde = (i - 1) * (K_TILDEBASE + 1) + K_TILDEBASE;

		// Point to bivGLWE(m / 2^{kappa_tilde * i})

		VecBiv* glwe_vec_ptr = ggsw->mat + row_i_ktilde * glwe_coef_number(params_glwe);

		// Computes the phase = m / 2^{kappa_tilde * i} + err
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec_ptr};
		glwe_secret_decrypt(module, phase_computed, sk_dft, &glwe_ct);

		// Computes the phase in Tn[X]
		biv_to_univ_rnx(params_glwe, phase_computed_univ_RnX, phase_computed);

		//! Computes by hand the phase = m / 2^{kappa_tilde * i}
		for (uint64_t p = 0; p < NBASE; p++)
		{
			phase_univ_RnX[p] = ldexp((double)m_univ[p], -(params_ggsw->kappa_tilde * i));
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
	delete_univ(phase_computed);
	delete_univ_rnx(phase_univ_RnX);
	delete_univ(m_univ);
	delete_univ_dft(m_univ_dft);
	delete_ggsw(ggsw);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_params(params_glwe);
	delete_ggsw_params(params_ggsw);
	pvda_delete_module_info(module);
}
