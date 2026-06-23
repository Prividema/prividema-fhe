#include "glwegadget_arithmetic.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"
#include "logger.h"
#include "maths_structures.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int glwegadget_half_prod(const MODULE* module, GLWECiphertext* result,
                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBiv* a)
{
	int status = -1;

	size_t nrows     = glwegadget_prep_ct->params->l_tilde;
	uint64_t nn      = glwegadget_prep_ct->params->params_glwe->nn;
	size_t ncols_in  = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);
	size_t ncols_out = glwe_params_n_limbs(result->params);

	GLWECiphertextDFT* glwe_dft = new_glwe_dft(result->params);
	CHECK_ALLOC(glwe_dft, "Allocation failed in half-product");

	CHECK_CALL(pvda_vmp_apply_dft(module, glwe_dft->vec, ncols_out, a, glwegadget_prep_ct->mat, nrows, ncols_in),
	           "vmp apply falied in half product");

	CHECK_CALL(glwe_dft_to_coef(module, result, glwe_dft),
	           "conversion from GLWE DFT to coefs failed in GLWEGadget half product");
	status = 0;
cleanup:
	delete_glwe_dft(glwe_dft);

	return status;
}

int glwegadget_half_prod_dft_to_dft(const MODULE* module, GLWECiphertextDFT* result_dft,
                                    const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBivDFT* a_dft)
{
	int status = -1;

	size_t nrows     = glwegadget_prep_ct->params->l_tilde;
	size_t ncols_in  = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);
	size_t ncols_out = glwe_params_n_limbs(result_dft->params);

	CHECK_CALL(pvda_vmp_apply_dft_to_dft(module, result_dft->vec, ncols_out, (double*)a_dft, nrows,
	                                     glwegadget_prep_ct->mat, nrows, ncols_in),
	           "vmp apply falied in half product");

	status = 0;
cleanup:
	return status;
}

int prepare_ksk(const MODULE* module, GLWEAutomorphismKSK* ksk, const GLWESecretKeyPrepared* new_key,
                const GLWESecretKeyPrepared* old_key)
{
	int status = -1;

	uint64_t k  = new_key->k;
	uint64_t nn = new_key->nn;

	GLWEGadgetCiphertext* glwegad_tmp = new_glwegadget(ksk->params);

	CHECK_ALLOC(glwegad_tmp, "GLWEGadget allocation failed in automorphism KSK preparation");

	ksk->automorphism_p = 0;

	for (int i = 0; i < k; ++i)
	{
		GLWEGadgetCiphertextPrep* gadget_ciph = ksk->enc_s[i];

		//GLWEGadget(sigma_p(sk_i))
		CHECK_CALL(
		    glwegadget_secret_encrypt(module, glwegad_tmp, new_key, glwe_prepared_sk_extract_poly_coefs(old_key, i)),
		    "GLWEGadget encryption failed in autmorphism KSK preparation");
		CHECK_CALL(glwegadget_prepare(module, gadget_ciph, glwegad_tmp),
		           "GLWEGadget preparation failed in automorphism KSK preparation");
	}

	status = 0;
cleanup:
	delete_glwegadget(glwegad_tmp);
	return status;
}
int prepare_automorphism_key(const MODULE* module, GLWEAutomorphismKSK* automorphism_ksk,
                             const GLWESecretKeyPrepared* glwe_key, int automorphism_p)
{
	int status = -1;

	if (!(automorphism_p & 1))
	{  //If autmorphism_p is even
		RAISE_ERROR("Cannot prepare autmorphism KSK for even p, operation is not well-defined");
	}

	uint64_t k  = glwe_key->k;
	uint64_t nn = glwe_key->nn;

	GLWEGadgetCiphertext* glwegad_tmp = new_glwegadget(automorphism_ksk->params);
	PolyUniv* auto_sk_tmp             = new_univ(automorphism_ksk->params->params_glwe);

	CHECK_ALLOC(glwegad_tmp, "GLWEGadget allocation failed in automorphism KSK preparation");
	CHECK_ALLOC(auto_sk_tmp, "Allocation failed in automorphism KSK preparation");

	automorphism_ksk->automorphism_p = automorphism_p;

	for (int i = 0; i < k; ++i)
	{
		GLWEGadgetCiphertextPrep* gadget_ciph = automorphism_ksk->enc_s[i];

		//sigma_p (sk_i)
		pvda_znx_automorphism(module, automorphism_p, auto_sk_tmp, glwe_prepared_sk_extract_poly_coefs(glwe_key, i));

		for (int i = 0; i < nn; ++i)
		{
			auto_sk_tmp[i] = -auto_sk_tmp[i];
		}

		//GLWEGadget(sigma_p(sk_i))
		CHECK_CALL(glwegadget_secret_encrypt(module, glwegad_tmp, glwe_key, auto_sk_tmp),
		           "GLWEGadget encryption failed in autmorphism KSK preparation");
		CHECK_CALL(glwegadget_prepare(module, gadget_ciph, glwegad_tmp),
		           "GLWEGadget preparation failed in automorphism KSK preparation");
	}

	status = 0;
cleanup:
	delete_glwegadget(glwegad_tmp);
	delete_univ(auto_sk_tmp);
	return status;
}

int glwegadget_automorphism(const MODULE* module, GLWECiphertext* result, const GLWEAutomorphismKSK* automorphism_ksk,
                            const GLWECiphertext* glwe)
{
	int status = -1;

	uint64_t nn            = result->params->nn;
	uint64_t k             = result->params->k;
	size_t nrows           = automorphism_ksk->params->l_tilde;
	uint64_t l_b_result    = glwe_params_l_b(result->params);
	int64_t automorphism_p = automorphism_ksk->automorphism_p;

	// This is the maximum internal precision of the result.
	// It is the maximum of the input b precision and the number of columns (GLWEGaget l_tilde precision) in the
	// key-switching key
	uint64_t biv_l = l_b_result > nrows ? l_b_result : nrows;
	if (k == 1)
	{
		PolyBiv* auto_tmp = new_biv_poly_custom_l(result->params, biv_l);
		CHECK_ALLOC(auto_tmp, "Allocation failed in automorphism");

		// auto_tmp = auto_p(a)
		PolyBiv a = glwe_extract_poly_view(glwe, 0);
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, &a);

		// result = halfProd(C_auto(-s), auto(a)) = -halfProd(C_auto(s), a)
		CHECK_CALL(glwegadget_half_prod(module, result, automorphism_ksk->enc_s[0], auto_tmp),
		           "half product in automorphism failed");

		// auto_tmp = auto_p(b)
		PolyBiv b = glwe_extract_poly_view(glwe, k);
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, &b);

		// result += auto_tmp ==> result = -halfProc(c_auto(s), auto(a)) + (0, auto(b))
		PolyBiv result_b = glwe_extract_poly_view(result, k);
		pvda_vec_znx_add(module, &result_b, &result_b, auto_tmp);

		status = 0;
	cleanup:
		delete_biv(auto_tmp);
	}
	else
	{
		PolyBiv* auto_tmp        = new_biv_poly_custom_l(result->params, biv_l);
		GLWECiphertext* glwe_tmp = new_glwe(result->params);

		CHECK_ALLOC_LABEL(auto_tmp, "Allocation failed in automorphism", cleanup2);
		CHECK_ALLOC_LABEL(glwe_tmp, "GLWE allocation failed in automorphism", cleanup2);

		// auto_tmp = auto_p(a_0)
		PolyBiv a_0 = glwe_extract_poly_view(glwe, 0);
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, &a_0);

		// result = halfProd(C_auto(-s_0), auto(a_0))
		CHECK_CALL_LABEL(glwegadget_half_prod(module, result, automorphism_ksk->enc_s[0], auto_tmp),
		                 "half product in automorphism failed", cleanup2);

		for (int i = 1; i < k; ++i)
		{
			// auto_tmp = auto_p(a_i)
			PolyBiv a_i = glwe_extract_poly_view(glwe, i);
			pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, &a_i);

			// result = halfProd(C_auto(-s_i), auto(a_i)) = -halfProd(C_auto(s_i), a_i)
			CHECK_CALL_LABEL(glwegadget_half_prod(module, glwe_tmp, automorphism_ksk->enc_s[i], auto_tmp),
			                 "half product in automorphism failed", cleanup2);

			add_glwe(module, result, result, glwe_tmp);
		}

		// auto_tmp = auto_p(b)
		PolyBiv b = glwe_extract_poly_view(glwe, k);
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, &b);

		// result += auto_tmp ==> result = -sum_i(halfProd(c_auto(s), auto(a_i))) + (0, auto(b))
		PolyBiv result_b = glwe_extract_poly_view(result, k);
		pvda_vec_znx_add(module, &result_b, &result_b, auto_tmp);
		status = 0;
	cleanup2:
		delete_biv(auto_tmp);
		delete_glwe(glwe_tmp);
	}

	return status;
}

int glwe_to_glwe_keyswitch(const MODULE* module, GLWECiphertext* result, const GLWEAutomorphismKSK* ksk,
                           const GLWECiphertext* glwe_ct)
{
	int status = -1;

	uint64_t k = result->params->k;

	if (k == 1)
	{
		PolyBiv a = glwe_extract_poly_view(glwe_ct, 0);
		//result = GLWE_k(k_new) \hp a
		//result = C_k(k_new) \hp a
		CHECK_CALL(glwegadget_half_prod(module, result, ksk->enc_s[0], &a), "half product in automorphism failed");

		// result = - (C_k(k_new) \hp a)
		negate_glwe(module, result, result);

		PolyBiv b        = glwe_extract_poly_view(glwe_ct, k);
		PolyBiv result_b = glwe_extract_poly_view(result, k);
		//result = (0, b) - result = (0, b) - (C_k(k_new) \hp a)
		pvda_vec_znx_add(module, &result_b, &result_b, &b);

		status = 0;
	cleanup:;
	}
	else
	{
		GLWECiphertext* glwe_tmp = new_glwe(result->params);

		CHECK_ALLOC_LABEL(glwe_tmp, "GLWE allocation failed in automorphism", cleanup2);

		PolyBiv a_0 = glwe_extract_poly_view(glwe_ct, 0);

		//result = GLWE_k(k_new[0]) \hp a_0
		//result = C_k(k_new[0]) \hp a_0
		CHECK_CALL_LABEL(glwegadget_half_prod(module, result, ksk->enc_s[0], &a_0),
		                 "half product in automorphism failed", cleanup2);

		for (int i = 1; i < k; ++i)
		{
			PolyBiv a_i = glwe_extract_poly_view(glwe_ct, i);

			//result += C_k(k_new[i]) \hp a_i
			//result = sum_{j=0}^{j=i} (C_k(k_new[j]) \hp a_j) (LOOP INVARIANT)
			CHECK_CALL_LABEL(glwegadget_half_prod(module, glwe_tmp, ksk->enc_s[i], &a_i),
			                 "half product in automorphism failed", cleanup2);

			add_glwe(module, result, result, glwe_tmp);
		}

		//result = -sum_{j=0}^{j=k-1} (C_k(k_new[j]) \hp a_j)
		negate_glwe(module, result, result);

		PolyBiv b        = glwe_extract_poly_view(glwe_ct, k);
		PolyBiv result_b = glwe_extract_poly_view(result, k);
		pvda_vec_znx_add(module, &result_b, &result_b, &b);
		status = 0;
	cleanup2:
		delete_glwe(glwe_tmp);
	}

	status = 0;
	return status;
}

int glwe_trace_expand(const MODULE* module, GLWECiphertext** results, int res_size, const GLWECiphertext* glwe_ct,
                      const GLWEAutomorphismKSKCollection* ksks)
{
	int status = -1;

	uint64_t nn = glwe_ct->params->nn;

	glwe_copy(results[0], glwe_ct);

	GLWECiphertext* tmp_glwe  = new_glwe(glwe_ct->params);
	GLWECiphertext* tmp_glwe2 = new_glwe(glwe_ct->params);
	CHECK_ALLOC(tmp_glwe, "Temp memory alloc in trace expansion failed");
	CHECK_ALLOC(tmp_glwe2, "Temp memory alloc in trace expansion failed");

	// Step 0:
	glwegadget_automorphism(module, tmp_glwe, glwegadget_ksk_collection_get_key(ksks, nn + 1), results[0]);

	add_glwe(module, tmp_glwe2, results[0], tmp_glwe);

	sub_glwe(module, tmp_glwe, results[0], tmp_glwe);
	PolyBiv tmp_flattened = glwe_flattened_biv(tmp_glwe);
	pvda_vec_znx_rotate(module, -1, &tmp_flattened, &tmp_flattened);
	normalize_glwe(module, results[1], tmp_glwe);
	normalize_glwe(module, results[0], tmp_glwe2);

	// Rest of the steps
	for (uint64_t p = 2; p < res_size; p *= 2)
	{
		int64_t auto_p = (int64_t)nn / p + 1;
		int64_t dist   = p;
		uint64_t b;
		for (b = 0; b < p && b + dist < res_size; ++b)
		{
			assert(b < res_size);
			GLWEAutomorphismKSK* ksk = glwegadget_ksk_collection_get_key(ksks, auto_p);
			CHECK_ALLOC(ksk, "KSK retrieval failed in trace expand");
			glwegadget_automorphism(module, tmp_glwe, ksk, results[b]);

			add_glwe(module, tmp_glwe2, results[b], tmp_glwe);

			sub_glwe(module, tmp_glwe, results[b], tmp_glwe);
			pvda_vec_znx_rotate(module, -p, &tmp_flattened, &tmp_flattened);
			normalize_glwe(module, results[b + dist], tmp_glwe);
			normalize_glwe(module, results[b], tmp_glwe2);
		}
		for (; b < p; ++b)
		{
			assert(b < res_size);
			GLWEAutomorphismKSK* ksk = glwegadget_ksk_collection_get_key(ksks, auto_p);
			CHECK_ALLOC(ksk, "KSK retrieval failed in trace expand");
			glwegadget_automorphism(module, tmp_glwe, ksk, results[b]);

			add_glwe(module, tmp_glwe, results[b], tmp_glwe);

			normalize_glwe(module, results[b], tmp_glwe);
		}
	}

	status = 0;
cleanup:
	delete_glwe(tmp_glwe);
	delete_glwe(tmp_glwe2);
	return status;
}

int packed_glwegadget_trace_expand(const MODULE* module, GLWEGadgetCiphertext** results, int res_size, int l_tilde,
                                   const GLWECiphertext* packed_glwegadget,
                                   const GLWEAutomorphismKSKCollection* auto_ksks)

{
	int status = -1;
	// TODO: add assertion/case to check if result has l_tilde less than l_tilde itself
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
			tmp->vec                                          = glwegadget_extract_bivglwe(results[res_num], prec_lvl);
			results_glwe[(prec_lvl - 1) * res_size + res_num] = tmp;
		}
	}

	CHECK_CALL(glwe_trace_expand(module, results_glwe, res_size * l_tilde, packed_glwegadget, auto_ksks),
	           "glwegadget_trace_expand failed in a GGSW trace expansion");

	status = 0;
cleanup:
	for (uint64_t i = 0; i < res_size; ++i)
		for (uint64_t j = 0; j < l_tilde; ++j) free(results_glwe[i * l_tilde + j]);

	return status;
}
int packed_glwegadget_trace_expand_prepared(const MODULE* module, GLWEGadgetCiphertextPrep** results, int res_size,
                                            const GLWECiphertext* packed_glwegadget,
                                            const GLWEAutomorphismKSKCollection* auto_ksks)
{
	int status = -1;

	GLWEGadgetCiphertext** gadgets = calloc(res_size, sizeof(GLWEGadgetCiphertext*));
	CHECK_ALLOC(gadgets, "unprepared gadget allocation failed in prepared glwegadget trace expansion");
	const GLWEGadgetParams* params_glwegad = results[0]->params;
	for (int r = 0; r < res_size; ++r)
	{
		gadgets[r] = new_glwegadget(params_glwegad);
		CHECK_ALLOC(gadgets[r], "GLWEGadget allocation failed in trace expansion");
	}

	CHECK_CALL(packed_glwegadget_trace_expand(module, gadgets, res_size, params_glwegad->l_tilde, packed_glwegadget,
	                                          auto_ksks),
	           "GLWEGadget trace expansion failed");

	for (int r = 0; r < res_size; ++r)
	{
		if (!results[r])
		{
			results[r] = new_glwegadget_prep(params_glwegad);
			CHECK_ALLOC(results[r], "Prepared GLWEGadget allocation failed");
		}
		CHECK_CALL(glwegadget_prepare(module, results[r], gadgets[r]),
		           "GLWEGadget preparation failed in trace expansion");
		delete_glwegadget(gadgets[r]);
	}

	free(gadgets);
	return 0;
cleanup:
	if (gadgets)
	{
		for (int r = 0; r < res_size; ++r)
		{
			delete_glwegadget(gadgets[r]);
		}
	}
	free(gadgets);
	return status;
}
