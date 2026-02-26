#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#include "core/ggsw/ggsw.h"
#include "rng.h"

#define NBASE            1024
#define KBASE            1
#define KAPPABASE        4
#define NLIMBSBASE       (KBASE + 1) * 1
#define LBASE            NLIMBSBASE / (KBASE + 1)
#define SIGMABASE        -(LBASE / 2 + 1) * KAPPABASE

#define K_TILDEBASE      1
#define KAPPA_TILDEBASE  4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1) * 1
#define L_TILDEBASE      NLIMBS_TILDEBASE / (K_TILDEBASE + 1)
#define SIGMA_TILDEBASE  -3


/**
 * @brief Tests ggsw_secret_encrpyt
 *
 */
Test(ggsw_secret_encrypt, works)
{
	// GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99%
	// chance
	double err_length = ldexp(1.0, -(LBASE / 2) * KAPPABASE) + 3 * sigma;
	cr_log_info("error length = %e", err_length);

	GLWECtParams* params_glwe_for_ggsw = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GGSWCtParams* params_ggsw =
	    new_ggsw_ct_params(params_glwe_for_ggsw, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
	MODULE* module           = new_module_info_p(NBASE);

	GGSWCiphertext* ct_ggsw  = new_ggsw(params_ggsw, NULL);
	GGSWSecretKeyDFT* sk_dft = new_uniform_ggsw_secret_key_dft(module, KBASE, 3);

	// Message uniformly drawn
	PolyUniv* msg_univ        = new_uniform_random_vec(NBASE, KAPPABASE * LBASE);
	PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
	vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, NBASE);

	// Computes a bivGGSW(msg)
	ggsw_secret_encrypt(module, params_ggsw, ct_ggsw, sk_dft, msg_univ);

	for (int64_t i = 1; i <= L_TILDEBASE; i++) {
		for (int64_t j = 0; j < K_TILDEBASE; j++) {
			// The pointer to bivGLWE(-m * sk_j / (2^kappa_tilde)^(i+1))
			VecBiv* ct_glwe = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			memcpy(ct_glwe, ct_ggsw->mat + ((i - 1) * (K_TILDEBASE + 1) + j) * glwe_coef_number(params_glwe_for_ggsw),
			       glwe_bytes(params_glwe_for_ggsw));

			// Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1)) + err
			PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			glwe_secret_demasking_ggsw_lib(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe);

			// Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1) + err in RnX
			PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);

			// Computes by hand the phase = -m * sk_j / 2^kappa_tilde*(i+1)
			PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			PolyUnivDFT* m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
			PolyUniv* m_skj_univ        = malloc(poly_univ_bytes(params_glwe_for_ggsw));

			// Computes DFT(msg * sk_j)
			mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);

			// Computes -msg * sk_j
			for (int64_t p = 0; p < NBASE; p++) {
				m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
			}
			vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1);

			// Computes -msg * sk_j / (2^kappa_tilde)^(i+1)
			for (int64_t p = 0; p < NBASE; p++) {
				phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -(params_ggsw->kappa_tilde * (i + 1)));
			}

			// Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase,
			// are only different by an error of approximation and a gaussian error
			for (int64_t p = 0; p < NBASE; p++) {
				double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
				double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] +
				                floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);
				double err_length = ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE);

				int cond =
				    (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

				cr_assert(
				    cond,
				    "Equality failed at p = %ld with : \n-(msg * sk_j)[%ld] / (2^kappa_tilde)^%ld = %lf and "
				    "phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \n-(msg * sk_j)[%ld] / "
				    "(2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf",
				    p, p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
				    phase_univ_RnX_computed[p], ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE),
				    p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
				    phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]),
				    ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE));
			}

			free(phase_univ_RnX_computed);
			free(m_skj_univ_dft);
			free(m_skj_univ);
			free(ct_glwe);
			free(phase);
			free(phase_univ_RnX);
		}
		// The pointer to bivGLWE(msg / (2^kappa_tilde)^(i+1))
		VecBiv* ct_glwe = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
		memcpy(ct_glwe,
		       ct_ggsw->mat + ((i - 1) * (K_TILDEBASE + 1) + K_TILDEBASE) * glwe_coef_number(params_glwe_for_ggsw),
		       glwe_bytes(params_glwe_for_ggsw));

		// Computes the phase = m/2^kappa_tilde + err
		PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
		glwe_secret_demasking_ggsw_lib(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe);

		// Computes the phase in Rn[X]
		PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
		biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);

		// Computes by hand the phase = msg / (2^kappa_tilde)^(i+1)
		PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));

		// Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase,
		// are only different by an error of approximation and a gaussian error
		for (int64_t p = 0; p < NBASE; p++) {
			double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
			double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] +
			                floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);

			int cond =
			    (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

			cr_assert(cond,
			          "Equality failed at p = %ld with : \nmsg[%ld] / (2^kappa_tilde)^%ld = %lf and "
			          "phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \nmsg[%ld] / (2^kappa_tilde)^%ld = "
			          "%lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf",
			          p, p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
			          phase_univ_RnX_computed[p], ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE),
			          p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
			          phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]),
			          ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE));
		}
		free(ct_glwe);
		free(phase);
		free(phase_univ_RnX);
		free(phase_univ_RnX_computed);
	}

	free(msg_univ);
	free(msg_univ_dft);
	delete_module_info_p(module);
	delete_ggsw(ct_ggsw);
	delete_ggsw_secret_key_dft(sk_dft);
	delete_glwe_ct_params(params_glwe_for_ggsw);
	delete_ggsw_ct_params(params_ggsw);
}

/**
 * @brief Construct a new Test object
 *
 */
Test(ggsw_secret_encrypt_dft, works)
{
	// GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
	double sigma = ldexp(1.0, -(LBASE / 2 + 1) * KAPPABASE);

	// The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99%
	// chance
	double err_length = ldexp(1.0, -(LBASE / 2) * KAPPABASE) + 3 * sigma;
	cr_log_info("error length = %e", err_length);

	GLWECtParams* params_glwe_for_ggsw = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
	GGSWCtParams* params_ggsw =
	    new_ggsw_ct_params(params_glwe_for_ggsw, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
	MODULE* module                 = new_module_info_p(NBASE);

	GGSWCiphertextDFT* ct_ggsw_dft = new_ggsw_dft(params_ggsw, NULL);
	GGSWSecretKeyDFT* sk_dft       = new_uniform_ggsw_secret_key_dft(module, KBASE, 3);

	// Message uniformly drawn
	PolyUniv* msg_univ        = new_uniform_random_vec(NBASE, KAPPABASE * LBASE);
	PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
	vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, NBASE);

	// Computes a bivGGSW(msg)
	ggsw_secret_encrypt_dft(module, params_ggsw, ct_ggsw_dft, sk_dft, msg_univ);

	for (int64_t i = 1; i <= L_TILDEBASE; i++) {
		for (int64_t j = 0; j < K_TILDEBASE; j++) {
			// The pointer to DFT(bivGLWE(-m * sk_j / (2^kappa_tilde)^(i+1)))
			VecBivDFT* ct_glwe_dft = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			memcpy(ct_glwe_dft,
			       ct_ggsw_dft->mat + ((i - 1) * (K_TILDEBASE + 1) + j) * glwe_coef_number(params_glwe_for_ggsw),
			       glwe_bytes(params_glwe_for_ggsw));

			// Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1)) + err
			PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			glwe_secret_demasking_ggsw_lib_dft(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe_dft);

			// Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1) + err in RnX
			PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);

			// Computes by hand the phase = -m * sk_j / 2^kappa_tilde*(i+1)
			PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
			PolyUnivDFT* m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
			PolyUniv* m_skj_univ        = malloc(poly_univ_bytes(params_glwe_for_ggsw));

			// Computes DFT(msg * sk_j)
			mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);

			// Computes -msg * sk_j
			for (int64_t p = 0; p < NBASE; p++) {
				m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
			}
			vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1);

			// Computes -msg * sk_j / (2^kappa_tilde)^(i+1)
			for (int64_t p = 0; p < NBASE; p++) {
				phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -(params_ggsw->kappa_tilde * (i + 1)));
			}

			// Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase,
			// are only different by an error of approximation and a gaussian error
			for (int64_t p = 0; p < NBASE; p++) {
				double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
				double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] +
				                floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);
				double err_length = ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE);

				int cond =
				    (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

				cr_assert(
				    cond,
				    "Equality failed at p = %ld with : \n-(msg * sk_j)[%ld] / (2^kappa_tilde)^%ld = %lf and "
				    "phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \n-(msg * sk_j)[%ld] / "
				    "(2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf",
				    p, p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
				    phase_univ_RnX_computed[p], ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE),
				    p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
				    phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]),
				    ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE));
			}

			free(phase_univ_RnX_computed);
			free(m_skj_univ_dft);
			free(m_skj_univ);
			free(ct_glwe_dft);
			free(phase);
			free(phase_univ_RnX);
		}
		// The pointer to bivGLWE(msg / (2^kappa_tilde)^(i+1))
		VecBivDFT* ct_glwe_dft = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
		memcpy(ct_glwe_dft,
		       ct_ggsw_dft->mat + ((i - 1) * (K_TILDEBASE + 1) + K_TILDEBASE) * glwe_coef_number(params_glwe_for_ggsw),
		       glwe_bytes(params_glwe_for_ggsw));

		// Computes the phase = m/2^kappa_tilde + err
		PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
		glwe_secret_demasking_ggsw_lib_dft(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe_dft);

		// Computes the phase in Rn[X]
		PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
		biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);

		// Computes by hand the phase = msg / (2^kappa_tilde)^(i+1)
		PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));

		// Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase,
		// are only different by an error of approximation and a gaussian error
		for (int64_t p = 0; p < NBASE; p++) {
			double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
			double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] +
			                floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);

			int cond =
			    (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

			cr_assert(cond,
			          "Equality failed at p = %ld with : \nmsg[%ld] / (2^kappa_tilde)^%ld = %lf and "
			          "phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \nmsg[%ld] / (2^kappa_tilde)^%ld = "
			          "%lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf",
			          p, p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
			          phase_univ_RnX_computed[p], ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE),
			          p, K_TILDEBASE, i + 1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p,
			          phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]),
			          ldexp(1.0, -(LBASE / 2) * KAPPABASE) + ldexp(1.0, -LBASE * KAPPABASE));
		}
		free(ct_glwe_dft);
		free(phase);
		free(phase_univ_RnX);
		free(phase_univ_RnX_computed);
	}

	free(msg_univ);
	free(msg_univ_dft);
	delete_module_info_p(module);
	delete_ggsw_dft(ct_ggsw_dft);
	delete_ggsw_secret_key_dft(sk_dft);
	delete_glwe_ct_params(params_glwe_for_ggsw);
	delete_ggsw_ct_params(params_ggsw);
}