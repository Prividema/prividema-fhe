#include "ggsw_ciphertext.h"

#include <string.h>

#include "glwe_ciphertext.h"
#include "logger.h"
#include "spqlios_alias.h"
#include "vec_znx_arithmetic_private.h"

//! GGSW Part (begin)

uint64_t ggsw_coef_number(const GGSWCtParams* params_ggsw)
{
	return params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe);
}

GGSWCiphertext* new_ggsw(const GGSWCtParams* params_ggsw)
{
	// The GGSW does not own the GGSWCtParams
	GGSWCiphertext* ggsw_ct = malloc(sizeof(GGSWCiphertext));
	if (log_is_null(ggsw_ct, "malloc in new_ggsw") < 0) return NULL;

	ggsw_ct->params = params_ggsw;

	// Initialize the GGSW ciphertext with 0s'
	ggsw_ct->mat = calloc(ggsw_coef_number(params_ggsw), sizeof(int64_t));
	if (log_is_null(ggsw_ct->mat, "calloc in new_ggsw") < 0)
	{
		free(ggsw_ct);
		return NULL;
	}

	return ggsw_ct;
}

void delete_ggsw(GGSWCiphertext* ggsw)
{
	free(ggsw->mat);
	free(ggsw);
}

VecBiv* ggsw_Sj_Yti(const GGSWCtParams* params_ggsw, MatBiv* ggsw_mat, int64_t j, int64_t i)
{
	// GLWE parameters
	uint64_t N       = params_ggsw->params_glwe->N;
	uint64_t k       = params_ggsw->params_glwe->k;
	uint64_t n_limbs = params_ggsw->params_glwe->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// GGSW parameters
	uint64_t k_tilde = params_ggsw->k_tilde;
	return ggsw_mat + (i - 1) * (k_tilde + 1) * n_limbs * N + j * n_limbs * N;
}

void normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw)
{
	// GGSW parameters
	const GGSWCtParams* params_ggsw = result->params;
	uint64_t k_tilde                = params_ggsw->k_tilde;
	uint64_t n_limbs_tilde          = params_ggsw->n_limbs_tilde;

	// GLWE parameters
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                      = params_glwe->N;
	uint64_t k                      = params_glwe->k;
	uint64_t n_limbs                = params_glwe->n_limbs;
	uint64_t kappa                  = params_glwe->kappa;
	uint64_t l                      = poly_biv_size(params_glwe);

	// Matrix parameters
	uint64_t nb_partial          = n_limbs_tilde / (k_tilde + 1);
	uint64_t nb_rows_per_partial = k_tilde + 1;

	for (uint64_t i = 1; i <= nb_partial; i++)
		for (uint64_t j = 0; j < nb_rows_per_partial; j++)
		{
			// The pointer to biGLWE(-m * sk_j * Y^i)
			// VecBiv* res_glwe = ggsw_Sj_Yti(params_ggsw, res->mat, j, i);
			VecBiv* res_glwe      = result->mat + (i - 1) * (k_tilde + 1) * n_limbs * N + j * n_limbs * N;
			const VecBiv* ct_glwe = ggsw_Sj_Yti(params_ggsw, ggsw->mat, j, i);

			// Normalize ct
			for (uint64_t t = 0; t < k + 1; t++)
				vec_znx_normalize_base2k_p(module, kappa, res_glwe + t * N, l, (k + 1) * N, ct_glwe + t * N, l,
				                           (k + 1) * N);
		}
}

void add_ggsw(GGSWCiphertext* res, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs)
{
	for (uint64_t t = 0; t < ggsw_coef_number(res->params); t++) 
		res->mat[t] = ggsw_lhs->mat[t] + ggsw_rhs->mat[t];
}

int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u_dft,
                    int do_normalization)
{
	// GGSW & GLWE params
	const GGSWCtParams* params_ggsw = result->params;
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t l       = poly_biv_size(params_glwe);
	int64_t mat_size = ggsw_size(params_ggsw);

	// The ciphertext in the DFT domain
	MatBivDFT* ct_dft = malloc(ggsw_bytes(params_ggsw));
	if (log_is_null(ct_dft, "malloc in const_mult_ggsw") < 0) return -1;

	vec_znx_dft_p(module, ct_dft, mat_size, ggsw->mat, mat_size, N);

	svp_apply_dft_p(module, ct_dft, mat_size, u_dft, ggsw->mat, mat_size, N);

	// Go back to Zn[X,Y]
	if (vec_znx_idft_p(module, result->mat, mat_size, ct_dft, mat_size) < 0)
		return log_perror("vec_znx_idft_p failed in const_mult_ggsw");

	// Normalization
	if (do_normalization)
		for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
			for (uint64_t j = 0; j < nb_rows_per_partial(params_ggsw); j++)
			{
				// The pointer to biGLWE(-m * sk_j / (2^{kappa_tilde}^i))
				VecBiv* ct_biv = ggsw_Sj_Yti(result->params, result->mat, j, i);
				for (uint64_t t = 0; t < k + 1; t++)
				{
					vec_znx_normalize_base2k_p(module, params_glwe->kappa, ct_biv + t * N, l, N * (k + 1),
					                           ct_biv + t * N, l, N * (k + 1));
				}
			}

	free(ct_dft);
	return 0;
}

//! GGSW DFT PART (begin)

uint64_t ggsw_coef_number_dft(const GGSWCtParams* params_ggsw)
{
	return (params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe)) / 2;
}

GGSWCiphertextDFT* new_ggsw_dft(const GGSWCtParams* params_ggsw)
{
	GGSWCiphertextDFT* ggsw_ct_dft = malloc(sizeof(GGSWCiphertextDFT));
	if (log_is_null(ggsw_ct_dft, "malloc in new_ggsw_dft") < 0) return NULL;

	ggsw_ct_dft->params = params_ggsw;

	// Initializes  the GGSW ciphertext with Os'
	ggsw_ct_dft->mat = calloc(2 * ggsw_coef_number_dft(params_ggsw), sizeof(double));
	if (log_is_null(ggsw_ct_dft->mat, "calloc in new_ggsw_dft") < 0)
	{
		free(ggsw_ct_dft);
		return NULL;
	}

	return ggsw_ct_dft;
}

void delete_ggsw_dft(GGSWCiphertextDFT* ggsw_dft)
{
	free(ggsw_dft->mat);
	free(ggsw_dft);
}

VecBivDFT* ggsw_Sj_Yti_dft(const GGSWCtParams* params_ggsw, MatBivDFT* ggsw_mat_dft, int64_t j, int64_t i)
{
	// GLWE parameters
	uint64_t N       = params_ggsw->params_glwe->N;
	uint64_t k       = params_ggsw->params_glwe->k;
	uint64_t n_limbs = params_ggsw->params_glwe->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// GGSW parameters
	uint64_t k_tilde = params_ggsw->k_tilde;

	return ggsw_mat_dft + (i - 1) * (k_tilde + 1) * n_limbs * N + j * n_limbs * N;
}

int normalize_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft)
{
	// GGSW parameters
	const GGSWCtParams* params_ggsw = ggsw_dft->params;
	uint64_t k_tilde                = params_ggsw->k_tilde;
	uint64_t n_limbs_tilde          = params_ggsw->n_limbs_tilde;

	// GLWE parameters
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                      = params_glwe->N;
	uint64_t k                      = params_glwe->k;
	uint64_t l                      = poly_biv_size(params_glwe);
	uint64_t n_limbs                = params_glwe->n_limbs;
	uint64_t kappa                  = params_glwe->kappa;

	// Matrix parameters
	uint64_t nb_partial          = n_limbs_tilde / (k_tilde + 1);
	uint64_t nb_rows_per_partial = k_tilde + 1;

	// The GGSW ciphertext's matrix out of DFT space
	MatBiv* ct_mat = malloc(ggsw_bytes(params_ggsw));
	if (log_is_null(ct_mat, "malloc in normalize_ggsw_dft") < 0) return -1;
	if (vec_znx_idft_p(module, ct_mat, ggsw_size(params_ggsw), ggsw_dft->mat, ggsw_size(params_ggsw)) < 0)
		return log_perror("vec_znx_idft_p failed in noramlize_ggsw_dft");

	for (uint64_t i = 1; i <= nb_partial; i++)
		for (uint64_t j = 0; j < nb_rows_per_partial; j++)
		{
			// The pointer to biGLWE(-m * sk_j * Y^i)
			VecBiv* ct_glwe = ggsw_Sj_Yti(params_ggsw, ct_mat, j, i);

			// Normalize ct
			// Normalize ct
			for (uint64_t t = 0; t < k + 1; t++)
				vec_znx_normalize_base2k_p(module, kappa, ct_glwe + t * N, l, (k + 1) * N, ct_glwe + t * N, l,
				                           (k + 1) * N);
		}

	// Computes the GGSW ciphertext's matrix in the DFT domain.
	vec_znx_dft_p(module, result_dft->mat, ggsw_size(params_ggsw), ct_mat, ggsw_size(params_ggsw), N);

	free(ct_mat);
	return 0;
}

void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft, const GGSWCiphertextDFT* ggsw_rhs_dft)
{
	for (uint64_t t = 0; t < ggsw_coef_number(result_dft->params); t++) 
		result_dft->mat[t] = ggsw_lhs_dft->mat[t] + ggsw_rhs_dft->mat[t];
}

int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft, int do_normalization)
{
	// GGSW & GLWE params
	const GGSWCtParams* params_ggsw = result_dft->params;
	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t l       = poly_biv_size(params_glwe);
	int64_t mat_size = ggsw_size(params_ggsw);

	// Temporary GGSW ciphertext
	MatBiv* tmp_ggsw_mat = malloc(ggsw_bytes(params_ggsw));
	if (log_is_null(tmp_ggsw_mat, "malloc in const_mult_ggsw_dft") < 0) return -1;

	// Computes tmp_ggsw = iDFT(ct_in_dft). Then computes :
	// ct_dft = DFT(u) * DFT(iDFT(ct_in_dft))) = DFT(u) * ct_in_dft
	if (vec_znx_idft_p(module, tmp_ggsw_mat, mat_size, ggsw_dft->mat, mat_size) < 0)
		return log_perror("vec_znx_idft_p failed in const_mult_ggsw_dft");
	svp_apply_dft_p(module, result_dft->mat, mat_size, u_dft, tmp_ggsw_mat, mat_size, N);

	// Normalization tmp_ggsw = u * iDFT(ct_dft_in)
	if (do_normalization)
	{
		// Computes res_dft out of DFT space
		if (vec_znx_idft_p(module, tmp_ggsw_mat, mat_size, result_dft->mat, mat_size) < 0)
			return log_perror("vec_znx_idft_p failed in const_mult_ggsw_dft");

		for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
			for (uint64_t j = 0; j < nb_rows_per_partial(params_ggsw); j++)
			{
				// The pointer to biGLWE(-m * sk_j * Y^i)
				VecBiv* ct_biv = ggsw_Sj_Yti(params_ggsw, tmp_ggsw_mat, j, i);
				for (uint64_t t = 0; t < k + 1; t++)
					vec_znx_normalize_base2k_p(module, params_glwe->kappa, ct_biv + t * N, l, N * (k + 1),
					                           ct_biv + t * N, l, N * (k + 1));
			}
		// Go back to DFT space
		vec_znx_dft_p(module, result_dft->mat, mat_size, tmp_ggsw_mat, mat_size, N);
	}

	free(tmp_ggsw_mat);
	return 0;
}

//! COMMON PART (begin)

uint64_t ggsw_size(const GGSWCtParams* params_ggsw) { return params_ggsw->n_limbs_tilde * params_ggsw->params_glwe->n_limbs; }

uint64_t ggsw_bytes(const GGSWCtParams* params_ggsw)
{
	int64_t N = params_ggsw->params_glwe->N;
	return ggsw_size(params_ggsw) * N * sizeof(int64_t);
}
