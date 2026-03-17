#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw.h"
#include "core/glwe/glwe.h"
#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_transform_key.h"
#include "rng.h"
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
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GLWESecretKey* sk        = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft = alloc_glwe_secret_key_dft(NBASE, KBASE);
	GGSWCiphertext* ggsw     = new_ggsw(params_ggsw);
	PolyUniv* m_univ         = malloc(poly_univ_bytes(params_glwe));
	PolyUnivDFT* m_univ_dft  = malloc(poly_univ_bytes(params_glwe));
	// Variables to compute the phase of each ggsw's row
	VecBiv* glwe_vec_computed            = calloc(glwe_coef_number(params_glwe), sizeof(int64_t));
	PolyBiv* phase_computed              = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* phase_computed_univ_RnX = calloc(NBASE, sizeof(int64_t));
	PolyUnivRnX* phase_univ_RnX          = calloc(NBASE, sizeof(int64_t));
	PolyUnivDFT* m_skj_univ_dft          = malloc(poly_univ_bytes(params_glwe));
	PolyUniv* m_skj_univ                 = malloc(poly_univ_bytes(params_glwe));

	// Draws uniformly in (Zn[X])^k the secret key
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Draws uniformly in Zn[X] the message
	uniform_random_pol_znx(m_univ, NBASE, KAPPABASE);

	// Computes the message in the DFT domain
	pvda_vec_znx_dft(module, m_univ_dft, 1, m_univ, 1, NBASE);

	// Computes a bivGGSW(m)
	ggsw_secret_encrypt(module, ggsw, sk_dft, m_univ);

	// Asserts the j-th row in the i-th PartialGGSW(m) in the ggsw is :
	// - a bivGLWE(-m * sk_j / 2^{kappa_tilde * i})), for j < k
	// - a bivGLWE(m / 2^{kappa_tilde * i}))        , for j = k
	for (uint64_t i = 1; i <= L_TILDEBASE; i++)
	{
		for (uint64_t j = 0; j < K_TILDEBASE; j++)
		{
			// Fills each changed variable with 0s'
			memset(phase_computed, 0, poly_biv_bytes(params_glwe));
			memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ_dft, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ, 0, poly_univ_bytes(params_glwe));
			memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

			// Copy bivGLWE(-m * sk_j / 2^{kappa_tilde * i}) in glwe_vec
			memcpy(glwe_vec_computed, ggsw_retrieve_bivglwe(params_ggsw, ggsw->mat, j, i), glwe_bytes(params_glwe));

			// Computes the phase = -m * sk_j / 2^{kappa_tilde * i}) + err
			GLWECiphertext glwe_ct = {params_glwe, glwe_vec_computed};
			glwe_secret_demasking(module, phase_computed, sk_dft, &glwe_ct);

			// Computes the phase = -m * sk_j / 2^{kappa_tilde * i} + err in RnX
			biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

			//! Computes by hand the phase = -m * sk_j / 2^{kappa_tilde*i}
			// Computes DFT(m * sk_j)
			mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, m_univ_dft, 1);

			// Computes -m * sk_j
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

		// Fills each changed variable with 0s'
		memset(phase_computed, 0, poly_biv_bytes(params_glwe));
		memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
		memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

		// The bivGGSW row correponding to bivGLWE(m / 2^{kappa_tilde * i})
		uint64_t row_i_ktilde = (i - 1) * (K_TILDEBASE + 1) + K_TILDEBASE;

		// Point to bivGLWE(m / 2^{kappa_tilde * i})
		memcpy(glwe_vec_computed, ggsw->mat + row_i_ktilde * glwe_coef_number(params_glwe), glwe_bytes(params_glwe));

		// Computes the phase = m / 2^{kappa_tilde * i} + err
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec_computed};
		glwe_secret_demasking(module, phase_computed, sk_dft, &glwe_ct);

		// Computes the phase in Tn[X]
		biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

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
	free(phase_computed_univ_RnX);
	free(m_skj_univ_dft);
	free(m_skj_univ);
	free(glwe_vec_computed);
	free(phase_computed);
	free(phase_univ_RnX);
	free(m_univ);
	free(m_univ_dft);
	delete_ggsw(ggsw);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_ct_params(params_glwe);
	delete_ggsw_ct_params(params_ggsw);
	pvda_delete_module_info(module);
}

/**
 * @brief Construct a new Test object
 *
 */
Test(ggsw_secret_encrypt_dft, works)
{
	// bivGLWE and bivGGSW parameters. This set of bivGLWE parameters is for bivGGSW ciphertext
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99%
	// chance
	double err_length = ldexp(1.0, -LBASE * KAPPABASE) + 3 * sigma;
	cr_log_info("error length = %e", err_length);

	// Parameters
	MODULE* module            = pvda_new_module_info(NBASE);
	GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

	// Variables
	GGSWCiphertextDFT* ggsw_dft = new_ggsw_dft(params_ggsw);
	GLWESecretKey* sk           = alloc_glwe_secret_key(NBASE, KBASE);
	GLWESecretKeyDFT* sk_dft    = alloc_glwe_secret_key_dft(NBASE, KBASE);
	PolyUniv* m_univ            = malloc(poly_univ_bytes(params_glwe));
	PolyUnivDFT* m_univ_dft     = malloc(poly_univ_bytes(params_glwe));

	// Variables to compute the phase of each ggsw's row
	VecBivDFT* glwe_vec_computed_dft     = calloc(glwe_coef_number(params_glwe), sizeof(int64_t));
	PolyBiv* phase_computed              = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	PolyUnivRnX* phase_computed_univ_RnX = calloc(NBASE, sizeof(int64_t));
	PolyUnivRnX* phase_univ_RnX          = calloc(NBASE, sizeof(int64_t));
	PolyUnivDFT* m_skj_univ_dft          = malloc(poly_univ_bytes(params_glwe));
	PolyUniv* m_skj_univ                 = malloc(poly_univ_bytes(params_glwe));

	// Draws uniformly in (Cm[X])^k the secret key in the DFT domain
	uniform_glwe_secret_key(module, sk, 3);
	transform_glwe_secret_key_not_dft_to_dft(module, sk_dft, sk);

	// Draws uniformly in Zn[X] the message
	uniform_random_pol_znx(m_univ, NBASE, KAPPABASE);

	// Computes the message in the DFT domain
	pvda_vec_znx_dft(module, m_univ_dft, 1, m_univ, 1, NBASE);

	// Computes a bivGGSW(m)
	ggsw_secret_encrypt_dft(module, ggsw_dft, sk_dft, m_univ);

	for (uint64_t i = 1; i <= L_TILDEBASE; i++)
	{
		for (uint64_t j = 0; j < K_TILDEBASE; j++)
		{
			// Fills each changed variable with 0s'
			memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
			memset(phase_computed, 0, poly_biv_bytes(params_glwe));
			memset(m_skj_univ_dft, 0, poly_univ_bytes(params_glwe));
			memset(m_skj_univ, 0, poly_univ_bytes(params_glwe));
			memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

			// The pointer to DFT(bivGLWE(-m * sk_j / 2^{kappa_tilde * i}))
			memcpy(glwe_vec_computed_dft, ggsw_retrieve_bivglwe_dft(params_ggsw, ggsw_dft->mat, j, i),
			       glwe_bytes(params_glwe));

			// Computes the phase = -m * sk_j / 2^{kappa_tilde * i}) + err
			GLWECiphertextDFT glwe_dft_ct = {params_glwe, glwe_vec_computed_dft};
			glwe_secret_demasking_dft(module, phase_computed, sk_dft, &glwe_dft_ct);

			// Computes the phase = -m * sk_j / 2^{kappa_tilde * i} + err in RnX
			biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

			//! Computes by hand the phase = -m * sk_j / 2^{kappa_tilde*i}
			// Computes DFT(m * sk_j)
			mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, m_univ_dft, 1);

			// Computes -m * sk_j
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

			// Assures that the difference between the phase = msg / (2^kappa_tilde)^i and the computed phase,
			// are only different by an error of approximation and a gaussian error
			for (uint64_t p = 0; p < NBASE; p++)
			{
				double diff = torus_distance(phase_univ_RnX[p], phase_computed_univ_RnX[p]);
				int cond    = diff < err_length;

				if (!cond) big_error_count++;
			}

			/// Proba that the number of error grater than 3sigma is greater or equal than 0.0027*N
			int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
			double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

			/// Asserts big_error_count <= 0.0027*N
			cr_assert(big_error_count <= max_fails,
			          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
			          "chance, that happens.",
			          max_fails, big_error_count, proba);
		}

		// Fills each changed variable with 0s'
		memset(phase_computed, 0, poly_biv_bytes(params_glwe));
		memset(phase_computed_univ_RnX, 0, poly_univ_bytes(params_glwe));
		memset(phase_univ_RnX, 0, poly_univ_bytes(params_glwe));

		// The bivGGSW row correponding to bivGLWE(m / 2^{kappa_tilde * i})
		uint64_t row_i_ktilde = ((i - 1) * (K_TILDEBASE + 1) + K_TILDEBASE);

		// The pointer to bivGLWE(m / 2^{kappa_tilde * i})
		memcpy(glwe_vec_computed_dft, ggsw_dft->mat + row_i_ktilde * glwe_coef_number(params_glwe),
		       glwe_bytes(params_glwe));

		// Computes the phase = m / 2^kappa_tilde + err
		GLWECiphertextDFT glwe_dft_ct = {params_glwe, glwe_vec_computed_dft};
		glwe_secret_demasking_dft(module, phase_computed, sk_dft, &glwe_dft_ct);

		// Computes the phase in Rn[X]
		biv_to_univ(params_glwe, phase_computed_univ_RnX, phase_computed);

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

		int max_fails = (int)(PROB_FACTOR * 0.0027 * NBASE);
		double proba  = binomial_tail(NBASE, 0.0027, PROB_FACTOR);

		/// Asserts big_error_count <= 0.0027*N
		cr_assert(big_error_count <= max_fails,
		          "The error should be greater than 3*sigma at most %ld times but got %ld times. There is a %lf "
		          "chance, that happens.",
		          max_fails, big_error_count, proba);
	}

	// Clean up
	free(phase_computed_univ_RnX);
	free(m_skj_univ_dft);
	free(m_skj_univ);
	free(glwe_vec_computed_dft);
	free(phase_computed);
	free(phase_univ_RnX);
	free(m_univ);
	free(m_univ_dft);
	delete_ggsw_dft(ggsw_dft);
	delete_glwe_secret_key_dft(sk_dft);
	delete_glwe_ct_params(params_glwe);
	delete_ggsw_ct_params(params_ggsw);
	pvda_delete_module_info(module);
}
