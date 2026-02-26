#include "ggsw_ciphertext.h"

#include "glwe_ciphertext.h"
#include "logger.h"
#include "spqlios_alias.h"
#include "vec_znx_arithmetic_private.h"

//! GGSW Part (begin)

uint64_t ggsw_coef_number(GGSWCtParams* params)
{
	return params->n_limbs_tilde * glwe_coef_number(params->params_glwe);
}

GGSWCiphertext* new_ggsw(GGSWCtParams* params, MatBiv* mat)
{
	GGSWCiphertext* ggsw_ct = malloc(sizeof(GGSWCiphertext));
	if (log_is_null(ggsw_ct, "malloc in new_ggsw") < 0) return NULL;

	ggsw_ct->params = params;

	if (mat == NULL) {
		ggsw_ct->mat = calloc(ggsw_coef_number(params), sizeof(int64_t));
		if (log_is_null(ggsw_ct->mat, "calloc in new_ggsw") < 0) {
			free(ggsw_ct);
			return NULL;
		}
	} else
		ggsw_ct->mat = mat;

	return ggsw_ct;
}

void delete_ggsw(GGSWCiphertext* ct)
{
	free(ct->mat);
	free(ct);
}

VecBiv* ggsw_Sj_Yti(GGSWCtParams* params_ggsw, MatBiv* ct_mat, int64_t j, int64_t i)
{
	// GLWE parameters
	uint64_t N       = params_ggsw->params_glwe->N;
	uint64_t k       = params_ggsw->params_glwe->k;
	uint64_t n_limbs = params_ggsw->params_glwe->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// GGSW parameters
	uint64_t k_tilde = params_ggsw->k_tilde;
	return ct_mat + (i - 1) * (k_tilde + 1) * n_limbs * N + j * n_limbs * N;
}

void normalize_ggsw(MODULE* module, GGSWCiphertext* res, GGSWCiphertext* ct)
{
	// GGSW parameters
	GGSWCtParams* params_ggsw = res->params;
	uint64_t k_tilde          = params_ggsw->k_tilde;
	uint64_t n_limbs_tilde    = params_ggsw->n_limbs_tilde;

	// GLWE parameters
	GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                = params_glwe->N;
	uint64_t k                = params_glwe->k;
	uint64_t n_limbs          = params_glwe->n_limbs;
	uint64_t kappa            = params_glwe->kappa;
	uint64_t l                = poly_biv_size(params_glwe);

	// Matrix parameters
	uint64_t nb_partial          = n_limbs_tilde / (k_tilde + 1);
	uint64_t nb_rows_per_partial = k_tilde + 1;

	for (int64_t i = 1; i <= nb_partial; i++)
		for (int64_t j = 0; j < nb_rows_per_partial; j++) {
			// The pointer to biGLWE(-m * sk_j * Y^i)
			// VecBiv* res_glwe = ggsw_Sj_Yti(params_ggsw, res->mat, j, i);
			VecBiv* res_glwe = res->mat + (i - 1) * (k_tilde + 1) * n_limbs * N + j * n_limbs * N;
			VecBiv* ct_glwe  = ggsw_Sj_Yti(params_ggsw, ct->mat, j, i);

			// Normalize ct
			for (int64_t t = 0; t < k + 1; t++)
				vec_znx_normalize_base2k_p(module, kappa, res_glwe + t * N, l, (k + 1) * N, ct_glwe + t * N, l,
				                           (k + 1) * N);
		}
}

void add_ggsw(GGSWCiphertext* res,  // result
              GGSWCiphertext* ct1,  // first operand
              GGSWCiphertext* ct2)  // second operand
{
	uint64_t nb_rows = res->params->n_limbs_tilde;
	uint64_t nb_cols = res->params->params_glwe->n_limbs;
	uint64_t N       = res->params->params_glwe->N;

	for (int64_t i = 0; i < nb_rows; i++)
		for (int64_t j = 0; j < nb_cols; j++)
			for (uint64_t k = 0; k < N; k++)
				res->mat[i * N * nb_cols + j * N + k] =
				    ct1->mat[i * N * nb_cols + j * N + k] + ct2->mat[i * N * nb_cols + j * N + k];
}

int const_mult_ggsw(MODULE* module, GGSWCiphertext* res, GGSWCiphertext* ct, PolyUnivDFT* u_dft, int do_normalization)
{
	// GGSW & GLWE params
	GGSWCtParams* params_ggsw = res->params;
	GLWECtParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N                = params_glwe->N;
	uint64_t k                = params_glwe->k;
	uint64_t l                = poly_biv_size(params_glwe);
	int64_t mat_size          = ggsw_size(params_ggsw);

	// The ciphertext in DFT space
	MatBivDFT* ct_dft = malloc(ggsw_bytes(params_ggsw));
	if (log_is_null(ct_dft, "malloc in const_mult_ggsw") < 0) return -1;

	vec_znx_dft_p(module, ct_dft, mat_size, ct->mat, mat_size, N);

	svp_apply_dft_p(module, ct_dft, mat_size, u_dft, ct->mat, mat_size, N);

	// Go back to Zn[X,Y]
	vec_znx_idft_p(module, res->mat, mat_size, ct_dft, mat_size);

	// Normalization
	if (do_normalization)
		for (int64_t i = 1; i <= nb_partials(params_ggsw); i++)
			for (int64_t j = 0; j < nb_rows_per_partial(params_ggsw); j++) {
				// The pointer to biGLWE(-m * sk_j / (2^{kappa_tilde}^i))
				VecBiv* ct_biv = ggsw_Sj_Yti(res->params, res->mat, j, i);
				for (int64_t t = 0; t < k + 1; t++) {
					vec_znx_normalize_base2k_p(module, params_glwe->kappa, ct_biv + t * N, l, N * (k + 1),
					                           ct_biv + t * N, l, N * (k + 1));
				}
			}

	free(ct_dft);
	return 0;
}

//! GGSW DFT PART (begin)

uint64_t ggsw_coef_number_dft(GGSWCtParams* params_ggsw)
{
	return (params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe)) / 2;
}

GGSWCiphertextDFT* new_ggsw_dft(GGSWCtParams* params, MatBivDFT* pmat)
{
	GGSWCiphertextDFT* ggsw_ct_dft = malloc(sizeof(GGSWCiphertextDFT));
	if (log_is_null(ggsw_ct_dft, "malloc in new_ggsw_dft") < 0) return NULL;

	ggsw_ct_dft->params = params;

	if (pmat == NULL) {
		ggsw_ct_dft->mat = calloc(2 * ggsw_coef_number_dft(params), sizeof(double));
		if (log_is_null(ggsw_ct_dft->mat, "calloc in new_ggsw_dft") < 0) {
			free(ggsw_ct_dft);
			return NULL;
		}
	} else
		ggsw_ct_dft->mat = pmat;

	return ggsw_ct_dft;
}

void delete_ggsw_dft(GGSWCiphertextDFT* res_dft)
{
	free(res_dft->mat);
	free(res_dft);
}

VecBivDFT* ggsw_Sj_Yti_dft(GGSWCtParams* params_ggsw, MatBivDFT* ct_dft, int64_t j, int64_t i)
{
	// GLWE parameters
	uint64_t N       = params_ggsw->params_glwe->N;
	uint64_t k       = params_ggsw->params_glwe->k;
	uint64_t n_limbs = params_ggsw->params_glwe->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// GGSW parameters
	uint64_t k_tilde = params_ggsw->k_tilde;

	return ct_dft + (i - 1) * (k_tilde + 1) * n_limbs * N + j * n_limbs * N;
}

int normalize_ggsw_dft(MODULE* module, GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct_dft)
{
	// GGSW parameters
	GGSWCtParams* params_ggsw = res_dft->params;
	uint64_t k_tilde          = params_ggsw->k_tilde;
	uint64_t n_limbs_tilde    = params_ggsw->n_limbs_tilde;

	// GLWE parameters
	GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                = params_glwe->N;
	uint64_t k                = params_glwe->k;
	uint64_t l                = poly_biv_size(params_glwe);
	uint64_t n_limbs          = params_glwe->n_limbs;
	uint64_t kappa            = params_glwe->kappa;

	// Matrix parameters
	uint64_t nb_partial          = n_limbs_tilde / (k_tilde + 1);
	uint64_t nb_rows_per_partial = k_tilde + 1;

	// The GGSW ciphertext's matrix out of DFT space
	MatBiv* ct_mat = malloc(ggsw_bytes(params_ggsw));
	if (log_is_null(ct_mat, "malloc in normalize_ggsw_dft") < 0) return -1;
	vec_znx_idft_p(module, ct_mat, ggsw_size(params_ggsw), ct_dft->mat, ggsw_size(params_ggsw));

	for (int64_t i = 1; i <= nb_partial; i++) 
		for (int64_t j = 0; j < nb_rows_per_partial; j++) {
			// The pointer to biGLWE(-m * sk_j * Y^i)
			VecBiv* ct_glwe = ggsw_Sj_Yti(params_ggsw, ct_mat, j, i);

			// Normalize ct
			// Normalize ct
			for (int64_t t = 0; t < k + 1; t++)
				vec_znx_normalize_base2k_p(module, kappa, ct_glwe + t * N, l, (k + 1) * N, ct_glwe + t * N, l,
				                           (k + 1) * N);
		}
	

	// Computes the GGSW ciphertext's matrix in DFT space.
	vec_znx_dft_p(module, res_dft->mat, ggsw_size(params_ggsw), ct_mat, ggsw_size(params_ggsw), N);

	free(ct_mat);
	return 0;
}

void add_ggsw_dft(GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct1_dft, GGSWCiphertextDFT* ct2_dft)
{
	uint64_t nb_rows = res_dft->params->n_limbs_tilde;
	uint64_t nb_cols = res_dft->params->params_glwe->n_limbs;
	uint64_t N       = res_dft->params->params_glwe->N;

	for (int64_t i = 0; i < nb_rows; i++)
		for (int64_t j = 0; j < nb_cols; j++)
			for (uint64_t k = 0; k < N; k++)
				res_dft->mat[i * N * nb_cols + j * N + k] =
				    ct1_dft->mat[i * N * nb_cols + j * N + k] + ct2_dft->mat[i * N * nb_cols + j * N + k];
}

int const_mult_ggsw_dft(MODULE* module, GGSWCiphertextDFT* res_dft, GGSWCiphertextDFT* ct_dft, PolyUnivDFT* u_dft,
                        int do_normalization)
{
	// GGSW & GLWE params
	GGSWCtParams* params_ggsw = res_dft->params;
	GLWECtParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N                = params_glwe->N;
	uint64_t k                = params_glwe->k;
	uint64_t l                = poly_biv_size(params_glwe);
	int64_t mat_size          = ggsw_size(params_ggsw);

	// Temporary GGSW ciphertext
	MatBiv* tmp_ggsw_mat = malloc(ggsw_bytes(params_ggsw));
	if (log_is_null(tmp_ggsw_mat, "malloc in const_mult_ggsw_dft") < 0) return -1;

	// Computes tmp_ggsw = iDFT(ct_in_dft). Then computes :
	// ct_dft = DFT(u) * DFT(iDFT(ct_in_dft))) = DFT(u) * ct_in_dft
	vec_znx_idft_p(module, tmp_ggsw_mat, mat_size, ct_dft->mat, mat_size);
	svp_apply_dft_p(module, res_dft->mat, mat_size, u_dft, tmp_ggsw_mat, mat_size, N);

	// Normalization tmp_ggsw = u * iDFT(ct_dft_in)
	if (do_normalization) {
		// Computes res_dft out of DFT space
		vec_znx_idft_p(module, tmp_ggsw_mat, mat_size, res_dft->mat, mat_size);

		for (int64_t i = 1; i <= nb_partials(params_ggsw); i++)
			for (int64_t j = 0; j < nb_rows_per_partial(params_ggsw); j++) {
				// The pointer to biGLWE(-m * sk_j * Y^i)
				VecBiv* ct_biv = ggsw_Sj_Yti(params_ggsw, tmp_ggsw_mat, j, i);
				for (int64_t t = 0; t < k + 1; t++)
					vec_znx_normalize_base2k_p(module, params_glwe->kappa, ct_biv + t * N, l, N * (k + 1),
					                           ct_biv + t * N, l, N * (k + 1));
			}
		// Go back to DFT space
		vec_znx_dft_p(module, res_dft->mat, mat_size, tmp_ggsw_mat, mat_size, N);
	}

	free(tmp_ggsw_mat);
	return 0;
}

//! COMMON PART (begin)

uint64_t ggsw_size(GGSWCtParams* params) { return params->n_limbs_tilde * params->params_glwe->n_limbs; }

uint64_t ggsw_bytes(GGSWCtParams* params)
{
	int64_t N = params->params_glwe->N;
	return ggsw_size(params) * N * sizeof(int64_t);
}
