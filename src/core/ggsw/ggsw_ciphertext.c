#include "ggsw_ciphertext.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "univariate_polynomial.h"
#include "utils.h"

// bivGGSW Part (begin)

GGSWCiphertext* new_ggsw(const GGSWParams* params_ggsw)
{
	// The bivGGSW does not own the GGSWParams
	GGSWCiphertext* ggsw = malloc(sizeof(GGSWCiphertext));
	CHECK_ALLOC(ggsw, "malloc in new_ggsw");

	ggsw->params = params_ggsw;

	// Initialize the bivGGSW ciphertext with 0s'
	ggsw->mat = calloc(ggsw_coef_number(params_ggsw), sizeof(int64_t));
	CHECK_ALLOC(ggsw->mat, "calloc in new_ggsw");

	return ggsw;
cleanup:

	free(ggsw);
	return NULL;
}

void delete_ggsw(GGSWCiphertext* ggsw)
{
	if (!ggsw) return;
	free(ggsw->mat);
	free(ggsw);
}

VecBiv* ggsw_retrieve_bivglwe(GGSWCiphertext* ggsw_ct, int64_t sk_idx, int64_t prec_lvl)
{
	assert(prec_lvl >= 1);

	// bivGLWE parameters
	const GLWEParams* params_glwe = ggsw_ct->params->params_glwe;

	// bivGGSW parameters
	uint64_t k_tilde = ggsw_ct->params->k_tilde;

	return ggsw_ct->mat + ((prec_lvl - 1) * (k_tilde + 1) + sk_idx) * glwe_coef_number(params_glwe);
}

int ggsw_secret_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWESecretKeyPrepared* sk_prep,
                        const PolyUniv* m_univ)
{
	int status = -1;

	const GGSWParams* params_ggsw = result->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;
	if (params_ggsw->k_tilde > params_ggsw->params_glwe->k)
		return log_perror("k_tilde should not be greater than k in ggsw_secret_encrypt");

	// bivGLWE parameters
	uint64_t nn      = params_glwe->nn;
	uint64_t k       = params_glwe->k;
	uint64_t k_tilde = params_ggsw->k_tilde;

	// Variables
	PolyUnivDFT* m_univ_dft     = new_univ_dft(module);   // DFT(msg)
	PolyUnivDFT* m_skj_univ_dft = new_univ_dft(module);   // DFT(msg * sk_j)
	PolyUniv* m_skj_univ        = new_univ(params_glwe);  // -msg * sk_j
	// Temp space for -m * sk * 2^{-kappa_tilde}
	PolyBiv* glwe_biv_msg = new_biv_poly(params_glwe);

	CHECK_ALLOC(m_univ_dft, "malloc failed in ggsw_secret_encrypt");
	CHECK_ALLOC(m_skj_univ_dft, "malloc failed in ggsw_secret_encrypt");
	CHECK_ALLOC(m_skj_univ, "malloc failed in ggsw_secret_encrypt");
	CHECK_ALLOC(glwe_biv_msg, "malloc failed in ggsw_secret_encrypt");

	// Computes DFT(msg)
	univ_coefs_to_dft(module, m_univ_dft, m_univ);

	for (uint64_t ij = 0; ij < ggsw_num_rows(params_ggsw); ++ij)
	{
		uint64_t j = (ij % (k + 1));
		uint64_t i = ij / (k + 1) + 1;
		{
			// For j < k, m_skj_univ is -m * skj
			// For j = k  we skip this and we get it from m_univ directly
			if (j < params_glwe->k)
			{
				// Computes DFT(msg * sk_j)
				mult_vec_znx_dft(module, m_skj_univ_dft, 1, glwe_prepared_sk_extract_poly_dft(sk_prep, j), 1,
				                 m_univ_dft, 1);

				// Computes -DFT(msg * sk_j)
				for (uint64_t p = 0; p < nn; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];

				// Invert the DFT to get -msg * sk_j
				CHECK_CALL(univ_dft_to_coefs(module, m_skj_univ, m_skj_univ_dft),
				           "vec_znx_idft_p failed in compute_phase_ij");
			}

			CHECK_CALL(univ_znx_to_biv(params_glwe, glwe_biv_msg, (k == j) ? m_univ : m_skj_univ,
			                           params_ggsw->kappa_tilde * i),
			           "univ_to_biv failed in compute_phase_ij");

			CHECK_CALL(add_biv_noise(module, params_glwe, glwe_biv_msg, glwe_biv_msg),
			           "Noise addition failed in GGSW encryption");
		}
		// Get the pointer for the result position
		VecBiv* glwe_vec       = ggsw_retrieve_bivglwe(result, j, i);
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec};

		//Compute: bivGLWE(glwe_biv_msg) into glwe_vec
		CHECK_CALL(glwe_secret_encrypt_phase(module, &glwe_ct, sk_prep, glwe_biv_msg),
		           "glwe_secret_masking_ggsw_lib failed in ggsw_secret_encrypt");
	}

	status = 0;

cleanup:
	free(glwe_biv_msg);
	delete_univ_dft(m_skj_univ_dft);
	delete_univ(m_skj_univ);
	delete_univ_dft(m_univ_dft);

	return status;
}
// bivGGSW DFT PART (begin)

GGSWCiphertextPrep* new_ggsw_prep(const GGSWParams* params_ggsw)
{
	GGSWCiphertextPrep* ggsw_mat_dft = malloc(sizeof(GGSWCiphertextPrep));
	CHECK_ALLOC(ggsw_mat_dft, "malloc in new_ggsw_dft");

	ggsw_mat_dft->params = params_ggsw;

	// Initializes  the bivGGSW ciphertext with Os'
	ggsw_mat_dft->mat = calloc(2 * ggsw_coef_number_dft(params_ggsw), sizeof(double));
	CHECK_ALLOC(ggsw_mat_dft->mat, "calloc in new_ggsw_dft");

	return ggsw_mat_dft;
cleanup:
	free(ggsw_mat_dft);
	return NULL;
}

void delete_ggsw_prep(GGSWCiphertextPrep* ggsw_dft)
{
	if (!ggsw_dft) return;
	free(ggsw_dft->mat);
	free(ggsw_dft);
}
int ggsw_prepare(const MODULE* module, GGSWCiphertextPrep* ggsw_prepared, const GGSWCiphertext* ggsw_ct)
{
	int status = -1;

	size_t nrows = ggsw_prepared->params->ciphertext_nb_limbs_tilde;
	size_t ncols = glwe_params_n_limbs(ggsw_prepared->params->params_glwe);

	CHECK_CALL(pvda_vmp_prepare_contiguous(module, ggsw_prepared->mat, ggsw_ct->mat, nrows, ncols),
	           "VMP prepare for GLWEGadget prepare failed");

	status = 0;
cleanup:
	return status;
}
