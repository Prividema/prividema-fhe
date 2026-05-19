#include "ggsw_arithmetic.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "bivariate_polynomial.h"
#include "ggsw_ciphertext.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwegadget_arithmetic.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "univariate_polynomial.h"
#include "utils.h"

// bivGGSW PART (begin)

int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw)
{
	int status                    = -1;
	const GGSWParams* params_ggsw = result->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;

	// Normalization of the bivGGSW ciphertext
	for (uint64_t limb_i = 0; limb_i < ggsw_num_rows(params_ggsw); ++limb_i)
	{
		uint64_t i = limb_i / (params_ggsw->k_tilde + 1) + 1;
		uint64_t j = (limb_i % (params_ggsw->k_tilde + 1));
		// The pointer to biGLWE(-m * sk_j * Y^i)
		VecBiv* result_glwe_vec       = ggsw_retrieve_bivglwe(result, j, i);
		GLWECiphertext glwe_normalize = {.params = params_glwe, .vec = result_glwe_vec};

		VecBiv* input_glwe_vec         = ggsw_retrieve_bivglwe(ggsw, j, i);
		GLWECiphertext input_normalize = {.params = ggsw->params->params_glwe, .vec = input_glwe_vec};

		normalize_glwe(module, &glwe_normalize, &input_normalize);
	}

	status = 0;

cleanup:

	return status;
}

void add_ggsw(const MODULE* module, GGSWCiphertext* res, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs)
{
	uint64_t nn = res->params->params_glwe->nn;
	pvda_vec_znx_add(module, res->mat, ggsw_total_n_glwe_limbs(res->params), nn, ggsw_lhs->mat,
	                 ggsw_total_n_glwe_limbs(ggsw_lhs->params), nn, ggsw_rhs->mat,
	                 ggsw_total_n_glwe_limbs(ggsw_rhs->params), nn);
}

int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u_dft)

{
	int status = -1;

	// Variables

	// bivGGSW & bivGLWE set of parameters
	const GGSWParams* params_ggsw = result->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;

	uint64_t nn                     = params_glwe->nn;
	int64_t mat_size                = ggsw_total_n_glwe_limbs(params_ggsw);
	GGSWCiphertextDFT* ggsw_tmp_dft = new_ggsw_dft(params_ggsw);

	CHECK_ALLOC(ggsw_tmp_dft, "alloc in const_mult_ggsw");

	pvda_svp_apply_dft(module, ggsw_tmp_dft->mat, mat_size, u_dft, ggsw->mat, mat_size, nn);

	// Go back to Zn[X,Y]
	CHECK_CALL(pvda_vec_znx_idft(module, result->mat, mat_size, ggsw_tmp_dft->mat, mat_size),
	           "vec_znx_idft_p failed in const_mult_ggsw");

	status = 0;

cleanup:
	delete_ggsw_dft(ggsw_tmp_dft);

	return status;
}

void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft,
                  const GGSWCiphertextDFT* ggsw_rhs_dft)
{
	for (uint64_t t = 0; t < ggsw_coef_number(result_dft->params); t++)
		result_dft->mat[t] = ggsw_lhs_dft->mat[t] + ggsw_rhs_dft->mat[t];
}

int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft)
{
	int status = -1;

	// Variables
	MatBiv* ggsw_mat = NULL;

	// bivGGSW & bivGLWEparams
	const GGSWParams* params_ggsw = result_dft->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;

	uint64_t nn       = params_glwe->nn;
	uint64_t mat_size = ggsw_total_n_glwe_limbs(params_ggsw);

	// Temporary bivGGSW ciphertext
	ggsw_mat = malloc(ggsw_bytes(params_ggsw));
	CHECK_ALLOC(ggsw_mat, "malloc in const_mult_ggsw_dft");

	// Computes ggsw_mat = iDFT(ggsw_mat_dft).
	CHECK_CALL(pvda_vec_znx_idft(module, ggsw_mat, mat_size, ggsw_dft->mat, mat_size),
	           "vec_znx_idft_p failed in const_mult_ggsw_dft");

	// Computes result_mat_dft = DFT(u) * DFT(iDFT(ggsw_mat_dft))) = DFT(u) * ggsw_mat_dft
	pvda_svp_apply_dft(module, result_dft->mat, mat_size, u_dft, ggsw_mat, mat_size, nn);

	status = 0;

cleanup:
	free(ggsw_mat);

	return status;
}

int ggsw_external_product(const MODULE* module,
                          GLWECiphertext* result,      // result
                          const GLWECiphertext* glwe,  // bivGLWE ciphertext
                          const GGSWCiphertext* ggsw   // bivGGSW ciphertext
)
{
	int status = -1;

	uint64_t nn = result->params->nn;
	// The bivGGSW ciphertext ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// The bivGLWE ciphertext glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// As the result of the vector-matrix product glwe * ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ggsw_num_rows(ggsw->params);
	uint64_t ncols = glwe_params_n_limbs(ggsw->params->params_glwe);

	MatBivDFT* ggsw_pmat  = NULL;  // Prepared bivGGSW ciphertext
	VecBivDFT* result_dft = NULL;  // ExternalProduct(glwe, ggsw)
	ggsw_pmat             = malloc(ggsw_bytes(ggsw->params));
	result_dft            = malloc(glwe_params_bytes(ggsw->params->params_glwe));

	CHECK_ALLOC(ggsw_pmat, "mat_dft's malloc failed in ggsw_external_product");

	CHECK_ALLOC(result_dft, "result's malloc failed in ggsw_external_product");

	CHECK_CALL(pvda_vmp_prepare_contiguous(module, ggsw_pmat, ggsw->mat, nrows, ncols),
	           "vmp_prepare_contiguous_p failed in ggsw_external_product");

	CHECK_CALL(pvda_vmp_apply_dft(module, result_dft, ncols, glwe->vec, glwe->params->ciphertext_nb_limbs, nn,
	                              ggsw_pmat, nrows, ncols),
	           "vmp_apply_dft_p failed in ggsw_external_product");

	CHECK_CALL(pvda_vec_znx_idft(module, result->vec, ncols, result_dft, ncols),
	           "vec_znx_idft_p failed in ggsw_external_product");

	status = 0;

cleanup:
	free(result_dft);
	free(ggsw_pmat);

	return 0;
}

int packed_glwegadget_trace_expand_ggsw(const MODULE* module, GGSWCiphertext** results, int res_size, int l_tilde,
                                        const GLWECiphertext* packed_glwegadget,
                                        const GLWEAutomorphismKSKCollection* auto_ksks,
                                        const GGSWCiphertext** sk_encryptions)
{
	assert(ggsw_params_l_tilde_a(results[0]->params) == l_tilde);
	assert(ggsw_params_l_tilde_b(results[0]->params) == l_tilde);
	int status = -1;
	GLWECiphertext* results_glwe[res_size * l_tilde];
	memset((uint8_t*)results_glwe, 0, sizeof(results_glwe));
	int64_t k = (int64_t)packed_glwegadget->params->k;

	/*
	 * Create dummy GLWECiphertext for the k'th GLWEs in each GGSW, that is, the ones that contain a
	 * m * 2^-jK.
	 * In other words, we are storing the results as GLWEGadgets inside the GGSWs by using only
	 * one every k GLWEs in a GGSW.
	 * That way, since we fill every k'th GLWE in a GGSW, to convert this "strided" GLWEGadget into a
	 * proper GGSW, it will suffice to generate the (-sk_i * m * 2^-jK) GLWEs that we are missing,
	 * which we can do by means of an external products with encryptions of -sk_i
	 */
	for (uint64_t prec_lvl = 1; prec_lvl <= l_tilde; ++prec_lvl)
	{
		for (uint64_t res_num = 0; res_num < res_size; ++res_num)
		{
			GLWECiphertext* tmp = malloc(sizeof(GLWECiphertext));
			CHECK_ALLOC(tmp, "Malloc failed in GGSW trace expansion");
			tmp->params                                       = results[res_num]->params->params_glwe;
			tmp->vec                                          = ggsw_retrieve_bivglwe(results[res_num], k, prec_lvl);
			results_glwe[(prec_lvl - 1) * res_size + res_num] = tmp;
		}
	}

	CHECK_CALL(glwe_trace_expand(module, results_glwe, res_size * l_tilde, packed_glwegadget, auto_ksks),
	           "glwegadget_trace_expand failed in a GGSW trace expansion");

	for (int res_num = 0; res_num < res_size; ++res_num)
	{
		for (uint64_t prec_lvl = 1; prec_lvl <= l_tilde; ++prec_lvl)
		{
			GLWECiphertext in = {results[res_num]->params->params_glwe,
			                     ggsw_retrieve_bivglwe(results[res_num], k, prec_lvl)};
			for (int i = 0; i < k; ++i)
			{
				GLWECiphertext res = {results[res_num]->params->params_glwe,
				                      ggsw_retrieve_bivglwe(results[res_num], i, prec_lvl)};
				CHECK_CALL(ggsw_external_product(module, &res, &in, sk_encryptions[i]),
				           "Relinearization external product failed in GGSW trace expansion");
				CHECK_CALL(normalize_glwe(module, &res, &res), "Normalization failed in GGSW trace expansion");
			}
		}
	}

	status = 0;
cleanup:
	for (uint64_t i = 0; i < res_size; ++i)
		for (uint64_t j = 0; j < l_tilde; ++j) free(results_glwe[i * l_tilde + j]);

	return status;
}
