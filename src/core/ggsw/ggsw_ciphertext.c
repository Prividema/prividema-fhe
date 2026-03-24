#include "ggsw_ciphertext.h"

#include <string.h>

#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "logger.h"
#include "spqlios_alias.h"
#include "utils.h"

//! bivGGSW Part (begin)

uint64_t ggsw_coef_number(const GGSWParams* params_ggsw)
{
	return params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe);
}

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

// TODO: add inline qualifier
VecBiv* ggsw_retrieve_bivglwe(GGSWCiphertext* ggsw_ct, int64_t j, int64_t i)
{
	// bivGLWE parameters
	const GLWEParams* params_glwe = ggsw_ct->params->params_glwe;

	// bivGGSW parameters
	uint64_t k_tilde = ggsw_ct->params->k_tilde;

	return ggsw_ct->mat + ((i - 1) * (k_tilde + 1) + j) * glwe_coef_number(params_glwe);
}

int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw)
{
	int status = -1;
	// TODO: assert input and output have equal params
	const GGSWParams* params_ggsw = result->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;
	uint64_t nn                   = params_glwe->nn;
	uint64_t k                    = params_glwe->k;
	uint64_t kappa                = params_glwe->kappa;
	uint64_t l                    = poly_biv_size(params_glwe);

	// Normalization of the bivGGSW ciphertext
	for (uint64_t i = 1; i <= ggsw_num_pggsw(params_ggsw); i++)
		for (uint64_t j = 0; j < ggsw_num_rows_per_pggsw(params_ggsw); j++)
		{
			// The pointer to biGLWE(-m * sk_j * Y^i)
			VecBiv* result_glwe_vec = ggsw_retrieve_bivglwe(result, j, i);

			// Normalize the k+1 bivGLWE's elements
			for (uint64_t t = 0; t < k + 1; t++)
				CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, result_glwe_vec + t * nn, l, (k + 1) * nn,
				                                         result_glwe_vec + t * nn, l, (k + 1) * nn),
				           "vec_znx_normalize_base2k_p failed in normalize_ggsw");
		}

	status = 0;

cleanup:

	return status;
}

void add_ggsw(GGSWCiphertext* res, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs)
{
	// TODO: move to spqlios
	for (uint64_t t = 0; t < ggsw_coef_number(res->params); t++) res->mat[t] = ggsw_lhs->mat[t] + ggsw_rhs->mat[t];
}

int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u_dft,
                    int do_normalization)
{
	int status = -1;

	// Variables
	MatBivDFT* ggsw_mat_dft = NULL;

	// bivGGSW & bivGLWE set of parameters
	const GGSWParams* params_ggsw = result->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N       = params_glwe->nn;
	int64_t mat_size = ggsw_size(params_ggsw);

	// The ciphertext in the DFT domain
	ggsw_mat_dft = malloc(ggsw_bytes(params_ggsw));
	CHECK_ALLOC(ggsw_mat_dft, "malloc in const_mult_ggsw");

	pvda_vec_znx_dft(module, ggsw_mat_dft, mat_size, ggsw->mat, mat_size, N);

	pvda_svp_apply_dft(module, ggsw_mat_dft, mat_size, u_dft, ggsw->mat, mat_size, N);

	// Go back to Zn[X,Y]
	CHECK_CALL(pvda_vec_znx_idft(module, result->mat, mat_size, ggsw_mat_dft, mat_size),
	           "vec_znx_idft_p failed in const_mult_ggsw");

	// Normalization
	if (do_normalization) normalize_ggsw(module, result, result);

	status = 0;

cleanup:
	free(ggsw_mat_dft);

	return status;
}

//! bivGGSW DFT PART (begin)

uint64_t ggsw_coef_number_dft(const GGSWParams* params_ggsw)
{
	return (params_ggsw->n_limbs_tilde * glwe_coef_number(params_ggsw->params_glwe)) / 2;
}

GGSWCiphertextDFT* new_ggsw_dft(const GGSWParams* params_ggsw)
{
	GGSWCiphertextDFT* ggsw_mat_dft = malloc(sizeof(GGSWCiphertextDFT));
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

void delete_ggsw_dft(GGSWCiphertextDFT* ggsw_dft)
{
	if (!ggsw_dft) return;
	free(ggsw_dft->mat);
	free(ggsw_dft);
}

VecBivDFT* ggsw_retrieve_bivglwe_dft(GGSWCiphertextDFT* ggsw_dft_ct, int64_t j, int64_t i)
{
	// bivGLWE parameters
	const GLWEParams* params_glwe = ggsw_dft_ct->params->params_glwe;

	// bivGGSW parameters
	uint64_t k_tilde = ggsw_dft_ct->params->k_tilde;

	return ggsw_dft_ct->mat + ((i - 1) * (k_tilde + 1) + j) * 2 * glwe_coef_number_dft(params_glwe);
}

int normalize_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft)
{
	//TODO: this function is untested and it is unclear that we actually want it implemented
	int status = -1;

	// Variables
	GGSWCiphertext* ggsw_ct = new_ggsw(ggsw_dft->params);
	CHECK_ALLOC(ggsw_ct, "Malloc failed in DFT GGSW normalization");

	// Computes the bivGGSW ciphertext out of the DFT domain
	CHECK_CALL(pvda_vec_znx_idft(module, ggsw_ct->mat, ggsw_size(result_dft->params), ggsw_dft->mat,
	                             ggsw_size(result_dft->params)),
	           "vec_znx_idft_p failed in noramlize_ggsw_dft");

	CHECK_CALL(normalize_ggsw(module, ggsw_ct, ggsw_ct), "Non-DFT GGSW normalization failed inside DFT normalization");

	// Computes the bivGGSW ciphertext's matrix in the DFT domain.
	pvda_vec_znx_dft(module, result_dft->mat, ggsw_size(result_dft->params), ggsw_ct->mat,
	                 ggsw_size(result_dft->params), result_dft->params->params_glwe->nn);

	status = 0;

cleanup:
	delete_ggsw(ggsw_ct);

	return status;
}

void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft,
                  const GGSWCiphertextDFT* ggsw_rhs_dft)
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
	const GGSWParams* params_ggsw = result_dft->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;

	uint64_t N       = params_glwe->nn;
	int64_t mat_size = ggsw_size(params_ggsw);

	// Temporary bivGGSW ciphertext
	ggsw_mat = malloc(ggsw_bytes(params_ggsw));
	CHECK_ALLOC(ggsw_mat, "malloc in const_mult_ggsw_dft");

	// Computes ggsw_mat = iDFT(ggsw_mat_dft). Then computes :
	CHECK_CALL(pvda_vec_znx_idft(module, ggsw_mat, mat_size, ggsw_dft->mat, mat_size),
	           "vec_znx_idft_p failed in const_mult_ggsw_dft");

	//TODO: check this, it is quite strange
	// Computes result_mat_dft = DFT(u) * DFT(iDFT(ggsw_mat_dft))) = DFT(u) * ggsw_mat_dft
	pvda_svp_apply_dft(module, result_dft->mat, mat_size, u_dft, ggsw_mat, mat_size, N);

	// Normalization of u * iDFT(ggsw_mat_dft)
	if (do_normalization) normalize_ggsw_dft(module, result_dft, result_dft);

	status = 0;

cleanup:
	free(ggsw_mat);

	return status;
}

//! COMMON PART (begin)

uint64_t ggsw_size(const GGSWParams* params_ggsw)
{
	return params_ggsw->n_limbs_tilde * params_ggsw->params_glwe->n_limbs;
}

uint64_t ggsw_bytes(const GGSWParams* params_ggsw)
{
	int64_t N = params_ggsw->params_glwe->nn;
	return ggsw_size(params_ggsw) * N * sizeof(int64_t);
}
