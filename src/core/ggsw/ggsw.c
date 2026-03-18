#include "ggsw.h"

#include <stdint.h>
#include <string.h>

#include "glwe.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwe_key.h"
#include "logger.h"
#include "math.h"
#include "rng.h"
#include "spqlios_alias.h"
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
	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t k_tilde = params_ggsw->k_tilde;

	// Variables
	PolyUnivDFT* m_univ_dft     = NULL;  // DFT(msg)
	PolyUnivDFT* m_skj_univ_dft = NULL;  // DFT(msg * sk_j)
	PolyUniv* m_skj_univ        = NULL;  // -msg * sk_j

	// compute_phase_ij requires some extra temp space
	// TODO: find out why it is being set to poly_univ_bytes even though it is a different data type
	double* tmp_sp1 = NULL;

	// Temp space for -m * sk * 2^{-kappa_tilde}
	PolyBiv* glwe_biv_msg = NULL;

	m_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_univ_dft, "malloc failed in ggsw_secret_encrypt");
	m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ_dft, "malloc failed in ggsw_secret_encrypt");
	m_skj_univ = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ, "malloc failed in ggsw_secret_encrypt");
	tmp_sp1 = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(tmp_sp1, "malloc failed in ggsw_secret_encrypt");
	glwe_biv_msg = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg, "malloc failed in ggsw_secret_encrypt");

	// Computes DFT(msg)
	pvda_vec_znx_dft(module, m_univ_dft, 1, m_univ, 1, N);

	for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
	{
		for (uint64_t j = 0; j < k_tilde + 1; j++)
		{
			{
				// For j < k, m_skj_univ is -m * skj
				// For j = k  we skip this and we get it from m_univ directly
				if (j < params_glwe->k)
				{
					// Computes DFT(msg * sk_j)
					mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, m_univ_dft, 1);

					// Computes -DFT(msg * sk_j)
					for (uint64_t p = 0; p < N; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];

					// Invert the DFT to get -msg * sk_j
					CHECK_CALL(pvda_vec_znx_idft(module, m_skj_univ, 1, m_skj_univ_dft, 1),
					           "vec_znx_idft_p failed in compute_phase_ij");
				}

				// Computes m_skj_univ / 2^{kappa_tilde*i}
				for (uint64_t p = 0; p < N; p++)
					tmp_sp1[p] =
					    ldexp((k == j) ? (double)m_univ[p] : (double)m_skj_univ[p], -params_ggsw->kappa_tilde * i);

				CHECK_CALL(add_normal_random_vec(tmp_sp1, N, tmp_sp1, 0.0, params_glwe->sigma),
				           "error addition failed in ggsw encryption");

				// Compute the base-2^kappa decomposition of tmp_sp1
				CHECK_CALL(univ_to_biv(params_glwe, glwe_biv_msg, tmp_sp1), "univ_to_biv failed in compute_phase_ij");
			}
			// Get the pointer for the result position
			VecBiv* glwe_vec       = ggsw_retrieve_bivglwe(params_ggsw, result->mat, j, i);
			GLWECiphertext glwe_ct = {params_glwe, glwe_vec};

			//Compute: bivGLWE(glwe_biv_msg) into glwe_vec
			CHECK_CALL(glwe_secret_masking(module, &glwe_ct, sk_dft, glwe_biv_msg),
			           "glwe_secret_masking_ggsw_lib failed in ggsw_secret_encrypt");
		}
	}

	status = 0;

cleanup:
	free(glwe_biv_msg);
	free(tmp_sp1);
	free(m_skj_univ_dft);
	free(m_skj_univ);
	free(m_univ_dft);

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
	uint64_t N = result->params->N;

	// The bivGGSW ciphertext ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// The bivGLWE ciphertext glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// As the result of the vector-matrix product glwe * ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ggsw->params->n_limbs_tilde;
	uint64_t ncols = ggsw->params->params_glwe->n_limbs;

	// Variables
	MatBivDFT* ggsw_pmat  = NULL;  // Prepared bivGGSW ciphertext
	VecBivDFT* result_dft = NULL;  // ExternalProduct(glwe, ggsw)

	ggsw_pmat = malloc(ggsw_bytes(ggsw->params));
	CHECK_ALLOC(ggsw_pmat, "mat_dft's malloc failed in ggsw_external_product");

	result_dft = malloc(glwe_bytes(ggsw->params->params_glwe));
	CHECK_ALLOC(result_dft, "result's malloc failed in ggsw_external_product");

	// Prepares bivGGSW ciphertext
	CHECK_CALL(pvda_vmp_prepare_contiguous(module, ggsw_pmat, ggsw->mat, nrows, ncols),
	           "vmp_prepare_contiguous_p failed in ggsw_external_product");

	// Computes ExternalProduct(glwe, ggsw)
	CHECK_CALL(pvda_vmp_apply_dft(module, result_dft, ncols, glwe->vec, nrows, N, ggsw_pmat, nrows, ncols),
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

//! bivGGSW IN DFT PART (begin)

int ggsw_secret_encrypt_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GLWESecretKeyDFT* sk_dft,
                            const PolyUniv* m_univ)
{
	int status = -1;

	const GGSWParams* params_ggsw = result_dft->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;
	if (params_ggsw->k_tilde > params_ggsw->params_glwe->k)
		return log_perror("k_tilde should not be greater than k in ggsw_secret_encrypt_dft");

	// bivGLWE parameters
	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t k_tilde = params_ggsw->k_tilde;

	// Variables
	PolyUnivDFT* m_univ_dft     = NULL;  // DFT(msg)
	PolyUnivDFT* m_skj_univ_dft = NULL;  // DFT(msg * sk_j)
	PolyUniv* m_skj_univ        = NULL;  // -m*sk_j

	double* msk_univ_RnX         = NULL;  // Temp space for -m * sk_j / 2^{...} (as univariate)
	PolyBiv* glwe_biv_msg        = NULL;  // Temp space for -m * sk_j / 2^{...} (as bivariate)
	PolyBivDFT* glwe_biv_msg_dft = NULL;  // We'll store DFT(-m * sk_j / 2^{kappa_tilde*i})

	m_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_univ_dft, "m_univ_dft's malloc failed");
	m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ_dft, "m_skj_univ_dft's malloc failed");
	m_skj_univ = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ, "m_skj_univ's alloc failed");
	msk_univ_RnX = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(msk_univ_RnX, "m_skj_univ's alloc failed");
	glwe_biv_msg = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg, "phase's malloc failed");
	glwe_biv_msg_dft = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg_dft, "phase_dft's malloc failed");

	// Computes DFT(m)
	pvda_vec_znx_dft(module, m_univ_dft, 1, m_univ, 1, N);

	for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
	{
		for (uint64_t j = 0; j < k_tilde + 1; j++)
		{
			// Computes: Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
			//           Dec_Kappa(m / 2^{kappa_tilde*i}) + err,         if j = k
			// The precision of the decomposition is l

			if (j < k)
			{
				// Computes DFT(msg * sk_j)
				mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, m_univ_dft, 1);

				// Computes -DFT(msg * sk_j)
				// TODO: study if accelerable using AVX, do we need it in spqlios?
				for (uint64_t p = 0; p < N; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];

				// Inverse DFT to retreive -msg * sk_j
				CHECK_CALL(pvda_vec_znx_idft(module, m_skj_univ, 1, m_skj_univ_dft, 1),
				           "vec_znx_idft_p failed in compute_phase_ij_dft");
			}
			for (uint64_t p = 0; p < N; p++)
				msk_univ_RnX[p] = ldexp((j == k) ? (double)m_univ[p] : m_skj_univ[p], -params_ggsw->kappa_tilde * i);

			// Add the error
			CHECK_CALL(add_normal_random_vec(msk_univ_RnX, N, msk_univ_RnX, 0.0, params_glwe->sigma),
			           "error addition failed in ggsw_dft encryption");

			// Convert the result to a bivariate (base-2k) polynomial
			CHECK_CALL(univ_to_biv(params_glwe, glwe_biv_msg, msk_univ_RnX),
			           "univ_to_biv failed in compute_phase_ij_dft");

			// Permorm DFT to get the result in the DFT domain
			pvda_vec_znx_dft(module, glwe_biv_msg_dft, poly_biv_size(params_glwe), glwe_biv_msg,
			                 poly_biv_size(params_glwe), N);

			// Result destination
			VecBivDFT* glwe_vec_dft = ggsw_retrieve_bivglwe_dft(params_ggsw, result_dft->mat, j, i);

			// Finally, mask/encrypt the above result to get a bivGLWE

			GLWECiphertextDFT glwe_dft_ct = {params_glwe, glwe_vec_dft};
			CHECK_CALL(glwe_secret_masking_dft(module, &glwe_dft_ct, sk_dft, glwe_biv_msg_dft),
			           "glwe_secret_masking_dft failed in ggsw_secret_encrypt_dft");
		}
	}

	status = 0;

cleanup:
	free(glwe_biv_msg_dft);
	free(glwe_biv_msg);
	free(msk_univ_RnX);
	free(m_skj_univ);
	free(m_skj_univ_dft);
	free(m_univ_dft);

	return 0;
}

int ggsw_external_product_dft(const MODULE* module,
                              GLWECiphertextDFT* result_dft,      // result
                              const GLWECiphertextDFT* glwe_dft,  // bivGLWE ciphertext
                              const GGSWCiphertextDFT* ggsw_dft   // bivGGSW ciphertext
)
{
	int status = -1;

	// Degree of chosen cyclotomic polynomial
	uint64_t N = result_dft->params->N;

	// The bivGLWE ciphertext glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// The bivGGSW ciphertext ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// As the result of the vector-matrix product glwe * ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ggsw_dft->params->n_limbs_tilde;
	uint64_t ncols = ggsw_dft->params->params_glwe->n_limbs;

	// Variables
	MatBiv* ggsw_mat     = NULL;
	MatBivDFT* ggsw_pmat = NULL;

	// Point to the bivGGSW ciphertext out of the DFT domain
	ggsw_mat = malloc(ggsw_bytes(ggsw_dft->params));
	CHECK_ALLOC(ggsw_mat, "mat's malloc failed in ggsw_external_product_dft");

	// Computes the bivGGSW ciphertext out of the DFT domain
	CHECK_CALL(pvda_vec_znx_idft(module, ggsw_mat, nrows * ncols, ggsw_dft->mat, nrows * ncols),
	           "vec_znx_idft_p failed in ggsw_external_product_dft");

	// Point to the bivGGSW ciphertext in the DFT domain
	ggsw_pmat = malloc(ggsw_bytes(ggsw_dft->params));
	CHECK_ALLOC(ggsw_pmat, "pmat's malloc failed in ggsw_external_product_dft");

	// Prepares the bivGGSW ciphertext in the DFT domain
	CHECK_CALL(pvda_vmp_prepare_contiguous(module, ggsw_pmat, ggsw_mat, nrows, ncols),
	           "vmp_prepare_contiguous_p failed in ggsw_external_product_dft");

	// Computes ExternalProduct(glwe, ggsw)
	CHECK_CALL(pvda_vmp_apply_dft_to_dft(module, result_dft->vec, ncols, glwe_dft->vec, nrows, ggsw_pmat, nrows, ncols),
	           "vmp_apply_dft_to_dft_p failed in ggsw_external_product_dft");

	status = 0;

cleanup:
	free(ggsw_pmat);
	free(ggsw_mat);

	return status;
}
