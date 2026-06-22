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
#include "maths_structures.h"
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
	uint64_t nn           = res->params->params_glwe->nn;
	PolyBiv res_flattened = ggsw_flattened_biv(res);
	PolyBiv lhs_flattened = ggsw_flattened_biv(ggsw_lhs);
	PolyBiv rhs_flattened = ggsw_flattened_biv(ggsw_rhs);
	pvda_vec_znx_add(module, &res_flattened, &lhs_flattened, &rhs_flattened);
}

int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* u_dft)

{
	int status = -1;

	// Variables

	// bivGGSW & bivGLWE set of parameters
	const GGSWParams* params_ggsw = result->params;
	const GLWEParams* params_glwe = params_ggsw->params_glwe;

	uint64_t nn                       = params_glwe->nn;
	int64_t mat_size                  = ggsw_total_n_glwe_limbs(params_ggsw);
	GGSWCiphertextPrep* ggsw_tmp_prep = new_ggsw_prep(params_ggsw);

	CHECK_ALLOC(ggsw_tmp_prep, "alloc in const_mult_ggsw");

	PolyBiv ggsw_flattened = ggsw_flattened_biv(ggsw);
	pvda_svp_apply_dft(module, ggsw_tmp_prep->mat, mat_size, u_dft, &ggsw_flattened);

	// Go back to Zn[X,Y]
	PolyBiv result_flattened = ggsw_flattened_biv(result);
	CHECK_CALL(pvda_vec_znx_idft(module, &result_flattened, ggsw_tmp_prep->mat, mat_size),
	           "vec_znx_idft_p failed in const_mult_ggsw");

	status = 0;

cleanup:
	delete_ggsw_prep(ggsw_tmp_prep);

	return status;
}

int ggsw_unprepared_external_product(const MODULE* module,
                                     GLWECiphertext* result,      // result
                                     const GLWECiphertext* glwe,  // bivGLWE ciphertext
                                     const GGSWCiphertext* ggsw   // bivGGSW ciphertext
)
{
	int status = -1;

	uint64_t nn = result->params->nn;

	uint64_t nrows     = ggsw_num_rows(ggsw->params);
	uint64_t ncols_in  = glwe_params_n_limbs(ggsw->params->params_glwe);
	uint64_t ncols_out = glwe_params_n_limbs(result->params);

	MatBivDFT* ggsw_pmat  = NULL;  // Prepared bivGGSW ciphertext
	VecBivDFT* result_dft = NULL;  // ExternalProduct(glwe, ggsw)
	ggsw_pmat             = malloc(ggsw_bytes(ggsw->params));
	result_dft            = malloc(glwe_params_bytes(ggsw->params->params_glwe));

	CHECK_ALLOC(ggsw_pmat, "mat_dft's malloc failed in ggsw_external_product");
	CHECK_ALLOC(result_dft, "result's malloc failed in ggsw_external_product");

	CHECK_CALL(pvda_vmp_prepare_contiguous(module, ggsw_pmat, ggsw->mat, nrows, ncols_in),
	           "vmp_prepare_contiguous_p failed in ggsw_external_product");

	PolyBiv glwe_flattened = glwe_flattened_biv(glwe);
	CHECK_CALL(pvda_vmp_apply_dft(module, result_dft, ncols_out, &glwe_flattened, ggsw_pmat, nrows, ncols_in),
	           "vmp_apply_dft_p failed in ggsw_external_product");

	PolyBiv result_flattened = glwe_flattened_biv(result);
	CHECK_CALL(pvda_vec_znx_idft(module, &result_flattened, result_dft, ncols_out),
	           "vec_znx_idft_p failed in ggsw_external_product");

	status = 0;

cleanup:
	free(result_dft);
	free(ggsw_pmat);

	return status;
}
int ggsw_external_product_to_dft(const MODULE* module, GLWECiphertextDFT* result, const GLWECiphertext* glwe,
                                 const GGSWCiphertextPrep* ggsw_prepared)
{
	int status = -1;

	uint64_t nn      = result->params->nn;
	uint64_t nrows   = ggsw_num_rows(ggsw_prepared->params);
	size_t ncols_in  = glwe_params_n_limbs(ggsw_prepared->params->params_glwe);
	size_t ncols_out = glwe_params_n_limbs(result->params);

	PolyBiv glwe_flattened = glwe_flattened_biv(glwe);
	CHECK_CALL(pvda_vmp_apply_dft(module, result->vec, ncols_out, &glwe_flattened, ggsw_prepared->mat, nrows, ncols_in),
	           "vmp_apply_dft_p failed in ggsw_external_product");

	status = 0;

cleanup:

	return status;
}
int ggsw_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                          const GGSWCiphertextPrep* ggsw_prepared)
{
	int status = -1;

	GLWECiphertextDFT* tmp_dft = new_glwe_dft(result->params);
	CHECK_ALLOC(tmp_dft, "Failed allocation in ggsw external product to non-dft");

	CHECK_CALL(ggsw_external_product_to_dft(module, tmp_dft, glwe, ggsw_prepared), "External product failed");
	CHECK_CALL(glwe_dft_to_coef(module, result, tmp_dft), "Failed iDFT in ggsw product");

	status = 0;
cleanup:
	delete_glwe_dft(tmp_dft);

	return status;
}

int packed_glwegadget_trace_expand_ggsw(const MODULE* module, GGSWCiphertext** results, int res_size, int l_tilde,
                                        const GLWECiphertext* packed_glwegadget,
                                        const GLWEAutomorphismKSKCollection* auto_ksks,
                                        const GGSWCiphertextPrep** sk_encryptions)
{
	assert(ggsw_params_l_tilde_a(results[0]->params) == l_tilde);
	assert(ggsw_params_l_tilde_b(results[0]->params) == l_tilde);
	int status = -1;
	GLWECiphertext* results_glwe[res_size * l_tilde];
	memset((uint8_t*)results_glwe, 0, sizeof(results_glwe));
	int64_t k = (int64_t)packed_glwegadget->params->k;

	/*
	 * Create dummy GLWECiphertext views for the k'th GLWE in each GGSW,
	 * that is, the ones that contain a m * 2^-jK.
	 * This way, glwe_trace_expand will fill the rows of the GGSWs that need to contain a m * 2^-jK.
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

	// Fills the k'th rows of GGSWs
	CHECK_CALL(glwe_trace_expand(module, results_glwe, res_size * l_tilde, packed_glwegadget, auto_ksks),
	           "glwegadget_trace_expand failed in a GGSW trace expansion");

	// Fills the rest of rows of the GGSWs, using GGSW(-s_i) encryptions to get
	// GLWE(-s_i * m * 2^-jK) from the GLWE(m * 2^-jK) we have
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

				// GGSW(-sk_i) HALFPROD GLWE(m * 2^-jK) = GLWE(-s_i * m * 2^-jK)
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

int packed_glwegadget_trace_expand_ggsw_prepared(const MODULE* module, GGSWCiphertextPrep** results, int res_size,
                                                 int l_tilde, const GLWECiphertext* packed_glwegadget,
                                                 const GLWEAutomorphismKSKCollection* auto_ksks,
                                                 const GGSWCiphertextPrep** sk_encryptions)
{
	int status = -1;

	GGSWCiphertext** ggsws = calloc(res_size, sizeof(GGSWCiphertextPrep*));
	CHECK_ALLOC(ggsws, "unprepared gadget allocation failed in prepared glwegadget trace expansion");
	const GGSWParams* params_ggsw = results[0]->params;
	for (int r = 0; r < res_size; ++r)
	{
		ggsws[r] = new_ggsw(params_ggsw);
		CHECK_ALLOC(ggsws[r], "GGSW allocation in trace expansion failed");
	}

	CHECK_CALL(packed_glwegadget_trace_expand_ggsw(module, ggsws, res_size, l_tilde, packed_glwegadget, auto_ksks,
	                                               sk_encryptions),
	           "GGSW trace expansion failed");

	for (int r = 0; r < res_size; ++r)
	{
		if (!results[r])
		{
			results[r] = new_ggsw_prep(params_ggsw);
			CHECK_ALLOC(results[r], "GGSW prepared allocation failed in trace expansion");
		}
		CHECK_CALL(ggsw_prepare(module, results[r], ggsws[r]), "GGSW preparation in trace expansion failed");
		delete_ggsw(ggsws[r]);
	}

	free(ggsws);
	return 0;
cleanup:
	if (ggsws)
	{
		for (int r = 0; r < res_size; ++r)
		{
			delete_ggsw(ggsws[r]);
		}
	}
	free(ggsws);
	return status;
}
