#include "ggsw.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "logger.h"
#include "math.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! bivGGSW PART (begin)

int ggsw_secret_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWESecretKeyDFT* sk_dft,
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
	// compute_phase_ij requires some extra temp space
	PolyUnivRnX* tmp_sp1 = new_univ_rnx(params_glwe);
	// Temp space for -m * sk * 2^{-kappa_tilde}
	PolyBiv* glwe_biv_msg = new_biv_poly(params_glwe);

	CHECK_ALLOC(m_univ_dft, "malloc failed in ggsw_secret_encrypt");
	CHECK_ALLOC(m_skj_univ_dft, "malloc failed in ggsw_secret_encrypt");
	CHECK_ALLOC(m_skj_univ, "malloc failed in ggsw_secret_encrypt");
	CHECK_ALLOC(tmp_sp1, "malloc failed in ggsw_secret_encrypt");
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
				mult_vec_znx_dft(module, m_skj_univ_dft, 1, glwe_sk_extract_poly_dft(sk_dft, j), 1, m_univ_dft, 1);

				// Computes -DFT(msg * sk_j)
				for (uint64_t p = 0; p < nn; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];

				// Invert the DFT to get -msg * sk_j
				CHECK_CALL(univ_dft_to_coefs(module, m_skj_univ, m_skj_univ_dft),
				           "vec_znx_idft_p failed in compute_phase_ij");
			}

			// Computes m_skj_univ / 2^{kappa_tilde*i}
			for (uint64_t p = 0; p < nn; p++)
				tmp_sp1[p] = ldexp((k == j) ? (double)m_univ[p] : (double)m_skj_univ[p], -params_ggsw->kappa_tilde * i);

			CHECK_CALL(add_normal_random_vec(tmp_sp1, nn, tmp_sp1, 0.0, params_glwe->sigma),
			           "error addition failed in ggsw encryption");

			// Compute the base-2^kappa decomposition of tmp_sp1
			CHECK_CALL(univ_rnx_to_biv(params_glwe, glwe_biv_msg, tmp_sp1), "univ_to_biv failed in compute_phase_ij");
		}
		// Get the pointer for the result position
		VecBiv* glwe_vec       = ggsw_retrieve_bivglwe(result, j, i);
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec};

		//Compute: bivGLWE(glwe_biv_msg) into glwe_vec
		CHECK_CALL(glwe_secret_encrypt_phase(module, &glwe_ct, sk_dft, glwe_biv_msg),
		           "glwe_secret_masking_ggsw_lib failed in ggsw_secret_encrypt");
	}

	status = 0;

cleanup:
	free(glwe_biv_msg);
	delete_univ_rnx(tmp_sp1);
	delete_univ_dft(m_skj_univ_dft);
	delete_univ(m_skj_univ);
	delete_univ_dft(m_univ_dft);

	return status;
}

int ggsw_external_product(const MODULE* module,
                          GLWECiphertext* result,      // result
                          const GLWECiphertext* glwe,  // bivGLWE ciphertext
                          const GGSWCiphertext* ggsw   // bivGGSW ciphertext
)
{
	int status = -1;

	// Degree of chosen cyclotomic polynomial
	uint64_t nn = result->params->nn;

	// The bivGGSW ciphertext ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// The bivGLWE ciphertext glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// As the result of the vector-matrix product glwe * ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ggsw_num_rows(ggsw->params);
	uint64_t ncols = glwe_params_n_limbs(ggsw->params->params_glwe);

	// Variables
	MatBivDFT* ggsw_pmat  = NULL;  // Prepared bivGGSW ciphertext
	VecBivDFT* result_dft = NULL;  // ExternalProduct(glwe, ggsw)
	                               //
	ggsw_pmat = malloc(ggsw_bytes(ggsw->params));
	CHECK_ALLOC(ggsw_pmat, "mat_dft's malloc failed in ggsw_external_product");

	result_dft = malloc(glwe_params_bytes(ggsw->params->params_glwe));
	CHECK_ALLOC(result_dft, "result's malloc failed in ggsw_external_product");

	// Prepares bivGGSW ciphertext
	CHECK_CALL(pvda_vmp_prepare_contiguous(module, ggsw_pmat, ggsw->mat, nrows, ncols),
	           "vmp_prepare_contiguous_p failed in ggsw_external_product");

	// Computes ExternalProduct(glwe, ggsw)
	CHECK_CALL(pvda_vmp_apply_dft(module, result_dft, ncols, glwe->vec, nrows, nn, ggsw_pmat, nrows, ncols),
	           "vmp_apply_dft_p failed in ggsw_external_product");

	// Computes the bivGGSW ciphertext out of the DFT domain
	CHECK_CALL(pvda_vec_znx_idft(module, result->vec, ncols, result_dft, ncols),
	           "vec_znx_idft_p failed in ggsw_external_product");

	status = 0;

cleanup:
	free(result_dft);
	free(ggsw_pmat);

	return 0;
}
