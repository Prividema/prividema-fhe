#include "ggsw_ciphertext.h"

#include <string.h>

#include "glwe_ciphertext.h"
#include "logger.h"
#include "utils.h"
#include "spqlios_alias.h"
#include "vec_znx_arithmetic_private.h"

//! bivGGSW Part (begin)

uint64_t ggsw_coef_number(const GGSWCtParams* params_ggsw)
{
	return params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe);
}

GGSWCiphertext* new_ggsw(const GGSWCtParams* params_ggsw)
{
	// The bivGGSW does not own the GGSWCtParams
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

// TODO: add inline qualifier
VecBiv* ggsw_retreive_bivglwe(const GGSWCtParams* params_ggsw, MatBiv* ggsw_mat, int64_t j, int64_t i)
{
	// bivGLWE parameters
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	// bivGGSW parameters
	uint64_t k_tilde = params_ggsw->k_tilde;

	return ggsw_mat + ((i - 1) * (k_tilde + 1) + j) * glwe_coef_number(params_glwe);
}

int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw)
{
	int status = -1;

	// bivGGSW parameters
	const GGSWCtParams* params_ggsw = result->params;
	uint64_t k_tilde                = params_ggsw->k_tilde;
	uint64_t n_limbs_tilde          = params_ggsw->n_limbs_tilde;

	// bivGLWE parameters
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                      = params_glwe->N;
	uint64_t k                      = params_glwe->k;
	uint64_t n_limbs                = params_glwe->n_limbs;
	uint64_t kappa                  = params_glwe->kappa;
	uint64_t l                      = poly_biv_size(params_glwe);

	// Matrix parameters
	uint64_t nb_partial          = n_limbs_tilde / (k_tilde + 1);
	uint64_t nb_rows_per_partial = k_tilde + 1;

	// Normalization of the bivGGSW ciphertext
	for (uint64_t i = 1; i <= nb_partial; i++)
		for (uint64_t j = 0; j < nb_rows_per_partial; j++)
		{
			// The pointer to biGLWE(-m * sk_j * Y^i)
			// VecBiv* res_glwe = ggsw_Sj_Yti(params_ggsw, res->mat, j, i);
			VecBiv* result_glwe_vec = ggsw_retreive_bivglwe(params_ggsw, result->mat, j, i);
			const VecBiv* ct_glwe = ggsw_retreive_bivglwe(params_ggsw, ggsw->mat, j, i);

			// Normalize the k+1 bivGLWE's elements
			for (uint64_t t = 0; t < k + 1; t++)
				CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, result_glwe_vec + t * N, l, (k + 1) * N, result_glwe_vec + t * N, l,
				                           (k + 1) * N), "vec_znx_normalize_base2k_p failed in normalize_ggsw");
		}

	status = 0;

cleanup:

	return status;
}

void add_ggsw(GGSWCiphertext* res, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs)
{
	for (uint64_t t = 0; t < ggsw_coef_number(res->params); t++) 
		res->mat[t] = ggsw_lhs->mat[t] + ggsw_rhs->mat[t];
}

int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u_dft,
                    int do_normalization)
{
	int status = -1;

	// Variables
	MatBivDFT* ggsw_mat_dft = NULL;

	// bivGGSW & bivGLWE set of parameters
	const GGSWCtParams* params_ggsw = result->params;
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t l       = poly_biv_size(params_glwe);
	int64_t mat_size = ggsw_size(params_ggsw);

	// The ciphertext in the DFT domain
	ggsw_mat_dft = malloc(ggsw_bytes(params_ggsw));
	CHECK_ALLOC(ggsw_mat_dft, "malloc in const_mult_ggsw");

	pvda_vec_znx_dft(module, ggsw_mat_dft, mat_size, ggsw->mat, mat_size, N);

	pvda_svp_apply_dft(module, ggsw_mat_dft, mat_size, u_dft, ggsw->mat, mat_size, N);

	// Go back to Zn[X,Y]
	CHECK_CALL(pvda_vec_znx_idft(module, result->mat, mat_size, ggsw_mat_dft, mat_size), "vec_znx_idft_p failed in const_mult_ggsw");

	// Normalization
	if (do_normalization)
		for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
			for (uint64_t j = 0; j < nb_rows_per_partial(params_ggsw); j++)
			{
				// The pointer to biGLWE(-m * sk_j / (2^{kappa_tilde}^i))
				VecBiv* glwe_vec = ggsw_retreive_bivglwe(result->params, result->mat, j, i);
				for (uint64_t t = 0; t < k + 1; t++)
					CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params_glwe->kappa, glwe_vec + t * N, l, N * (k + 1),
					                           glwe_vec + t * N, l, N * (k + 1)), "vec_znx_normalize_base2k_p failed in const_mult_ggsw");
			}

	status = 0;

cleanup:
	free(ggsw_mat_dft);
	
	return status;
}

//! bivGGSW DFT PART (begin)

uint64_t ggsw_coef_number_dft(const GGSWCtParams* params_ggsw)
{
	return (params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe)) / 2;
}

GGSWCiphertextDFT* new_ggsw_dft(const GGSWCtParams* params_ggsw)
{
	GGSWCiphertextDFT* ggsw_mat_dft = malloc(sizeof(GGSWCiphertextDFT));
  CHECK_ALLOC(ggsw_mat_dft,"malloc in new_ggsw_dft");

	ggsw_mat_dft->params = params_ggsw;

	// Initializes  the bivGGSW ciphertext with Os'
	ggsw_mat_dft->mat = calloc(2 * ggsw_coef_number_dft(params_ggsw), sizeof(double));
  CHECK_ALLOC(ggsw_mat_dft->mat, "calloc in new_ggsw_dft");

	return ggsw_mat_dft;
cleanup:
  free(ggsw_mat_dft);
  return NULL;
}

void delete_ggsw_dft(GGSWCiphertextDFT* ggsw_dft)
{
  if (!ggsw_dft) return;
	free(ggsw_dft->mat);
	free(ggsw_dft);
}

VecBivDFT* ggsw_retreive_bivglwe_dft(const GGSWCtParams* params_ggsw, MatBivDFT* ggsw_mat_dft, int64_t j, int64_t i)
{
	// bivGLWE parameters
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	// bivGGSW parameters
	uint64_t k_tilde = params_ggsw->k_tilde;

	return ggsw_mat_dft + ((i - 1) * (k_tilde + 1) + j) * 2 * glwe_coef_number_dft(params_glwe);
}

int normalize_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft)
{
	int status = -1;

	// Variables 
	MatBiv* ggsw_mat = NULL;

	// bivGGSW parameters
	const GGSWCtParams* params_ggsw = ggsw_dft->params;
	uint64_t k_tilde                = params_ggsw->k_tilde;
	uint64_t n_limbs_tilde          = params_ggsw->n_limbs_tilde;

	// bivGLWE parameters
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                      = params_glwe->N;
	uint64_t k                      = params_glwe->k;
	uint64_t l                      = poly_biv_size(params_glwe);
	uint64_t n_limbs                = params_glwe->n_limbs;
	uint64_t kappa                  = params_glwe->kappa;

	// Matrix parameters
	uint64_t nb_partial          = n_limbs_tilde / (k_tilde + 1);
	uint64_t nb_rows_per_partial = k_tilde + 1;

	// Points to the bivGGSW ciphertext out of the DFT domain
	ggsw_mat = malloc(ggsw_bytes(params_ggsw));
	CHECK_ALLOC(ggsw_mat, "malloc in normalize_ggsw_dft");

	// Computes the bivGGSW ciphertext out of the DFT domain
	CHECK_CALL(pvda_vec_znx_idft(module, ggsw_mat, ggsw_size(params_ggsw), ggsw_dft->mat, ggsw_size(params_ggsw)), 
		"vec_znx_idft_p failed in noramlize_ggsw_dft");

	// Normalization of the bivGGSW ciphertext
	for (uint64_t i = 1; i <= nb_partial; i++)
		for (uint64_t j = 0; j < nb_rows_per_partial; j++)
		{
			// The pointer to biGLWE(-m * sk_j * Y^i)
			VecBiv* glwe_vec = ggsw_retreive_bivglwe(params_ggsw, ggsw_mat, j, i);

			// Normalize ct
			for (uint64_t t = 0; t < k + 1; t++)
				CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, glwe_vec + t * N, l, (k + 1) * N, glwe_vec + t * N, l,
				                           (k + 1) * N), "vec_normalize_base2k_p failed in normalize_ggsw_dft");
		}

	// Computes the bivGGSW ciphertext's matrix in the DFT domain.
	pvda_vec_znx_dft(module, result_dft->mat, ggsw_size(params_ggsw), ggsw_mat, ggsw_size(params_ggsw), N);

	status = 0;

cleanup:
	free(ggsw_mat);
	
	return status;
}

void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft, const GGSWCiphertextDFT* ggsw_rhs_dft)
{
	for (uint64_t t = 0; t < ggsw_coef_number(result_dft->params); t++) 
		result_dft->mat[t] = ggsw_lhs_dft->mat[t] + ggsw_rhs_dft->mat[t];
}

int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft, int do_normalization)
{
	int status = -1;

	// Variables 
	MatBiv* ggsw_mat = NULL;

	// bivGGSW & bivGLWEparams
	const GGSWCtParams* params_ggsw = result_dft->params;
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t l       = poly_biv_size(params_glwe);
	int64_t mat_size = ggsw_size(params_ggsw);

	// Temporary bivGGSW ciphertext
	ggsw_mat = malloc(ggsw_bytes(params_ggsw));
	CHECK_ALLOC(ggsw_mat, "malloc in const_mult_ggsw_dft");

	// Computes ggsw_mat = iDFT(ggsw_mat_dft). Then computes :
	CHECK_CALL(pvda_vec_znx_idft(module, ggsw_mat, mat_size, ggsw_dft->mat, mat_size), "vec_znx_idft_p failed in const_mult_ggsw_dft");
	
	// Computes result_mat_dft = DFT(u) * DFT(iDFT(ggsw_mat_dft))) = DFT(u) * ggsw_mat_dft
	pvda_svp_apply_dft(module, result_dft->mat, mat_size, u_dft, ggsw_mat, mat_size, N);

	// Normalization of u * iDFT(ggsw_mat_dft)
	if (do_normalization)
	{
		// Computes ggsw_mat = u * iDFT(ggsw_mat_dft)
		CHECK_CALL(pvda_vec_znx_idft(module, ggsw_mat, mat_size, result_dft->mat, mat_size), "vec_znx_idft_p failed in const_mult_ggsw_dft");

		for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
			for (uint64_t j = 0; j < nb_rows_per_partial(params_ggsw); j++)
			{
				// The pointer to biGLWE(-m * sk_j * Y^i)
				VecBiv* glwe_vec = ggsw_retreive_bivglwe(params_ggsw, ggsw_mat, j, i);
				// Normalize the k+1 bivGLWE's elements
				for (uint64_t t = 0; t < k + 1; t++)
					CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params_glwe->kappa, glwe_vec + t * N, l, N * (k + 1),
					                           glwe_vec + t * N, l, N * (k + 1)), "vec_normalize_base2k_p failed in const_mult_ggsw_dft");
			}
		// Go back to the DFT domain
		pvda_vec_znx_dft(module, result_dft->mat, mat_size, ggsw_mat, mat_size, N);
	}

	status = 0;

cleanup:
	free(ggsw_mat);
	
	return status;
}

//! COMMON PART (begin)

uint64_t ggsw_size(const GGSWCtParams* params_ggsw) { return params_ggsw->n_limbs_tilde * params_ggsw->params_glwe->n_limbs; }

uint64_t ggsw_bytes(const GGSWCtParams* params_ggsw)
{
	int64_t N = params_ggsw->params_glwe->N;
	return ggsw_size(params_ggsw) * N * sizeof(int64_t);
}
