#include "glwegadget_arithmetic.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_arithmetic.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwegadget_ciphertext.h"
#include "logger.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

int glwegadget_half_prod(const MODULE* module, GLWECiphertext* result,
                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBiv* a)
{
	int status = -1;

	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	uint64_t nn  = glwegadget_prep_ct->params->params_glwe->nn;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	GLWECiphertextDFT* glwe_dft = new_glwe_dft(result->params);

	CHECK_CALL(pvda_vmp_apply_dft(module, glwe_dft->vec, ncols, a, nrows, nn, glwegadget_prep_ct->mat, nrows, ncols),
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

	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	CHECK_CALL(pvda_vmp_apply_dft_to_dft(module, result_dft->vec, ncols, (double*)a_dft, nrows, glwegadget_prep_ct->mat,
	                                     nrows, ncols),
	           "vmp apply falied in half product");

	status = 0;
cleanup:
	return status;
}

int prepare_automorphism_key(const MODULE* module, GLWEAutomorphismKSK* automorphism_ksk,
                             const GLWESecretKeyPrepared* glwe_key, int automorphism_p)
{
	int status = -1;

	uint64_t k  = glwe_key->k;
	uint64_t nn = glwe_key->nn;

	GLWEGadgetCiphertext* glwegad_tmp = new_glwegadget(automorphism_ksk->params);
	CHECK_ALLOC(glwegad_tmp, "GLWEGadget allocation failed in automorphism KSK preparation");
	PolyUniv* auto_sk_tmp = new_univ(automorphism_ksk->params->params_glwe);

	automorphism_ksk->automorphism_p = automorphism_p;

	for (int i = 0; i < k; ++i)
	{
		GLWEGadgetCiphertextPrep* gadget_ciph = automorphism_ksk->enc_s[i];

		//sigma_p (sk_i)
		pvda_vec_znx_automorphism(module, automorphism_p, auto_sk_tmp, 1, nn,
		                          glwe_prepared_sk_extract_poly_coefs(glwe_key, i), 1, nn);

		//GLWEGadget(sigma_p(sk_i))
		glwegadget_secret_encrypt(module, glwegad_tmp, glwe_key, auto_sk_tmp);
		glwegadget_prepare(module, gadget_ciph, glwegad_tmp);
	}

	status = 0;
cleanup:
	delete_glwegadget(glwegad_tmp);
	return status;
}

int glwegadget_automorphism(const MODULE* module, GLWECiphertext* result, const GLWEAutomorphismKSK* automorphism_ksk,
                            const GLWECiphertext* glwe, int automorphism_p)
{
	int status = -1;

	if (!(automorphism_p & 1))  // if automorphism_p  is even
	{
		log_message(LOG_WARN,
		            "Trying to perform an automorphism with even k. Program will continue with unpredictable results");
	}

	if (automorphism_p != automorphism_ksk->automorphism_p)
	{
		log_message(LOG_WARN,
		            "Trying to perform an automorphism with a KSK with mismatched k. Results might be unexpected");
	}

	uint64_t nn         = result->params->nn;
	uint64_t k          = result->params->k;
	size_t nrows        = automorphism_ksk->params->l_tilde;
	uint64_t l_b_result = glwe_params_l_b(result->params);

	// This is the maximum internal precision of the result.
	// It is the maximum of the input b precision and the number of columns (GLWEGaget l_tilde precision) in the
	// key-switching key
	uint64_t biv_l = l_b_result > nrows ? l_b_result : nrows;
	if (k == 1)
	{
		PolyBiv* auto_tmp = new_biv_poly_custom_l(result->params, biv_l);

		// auto_tmp = auto_p(a)
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, biv_l, nn, glwe->vec, glwe_params_l_a(glwe->params),
		                          2 * nn);

		// result = halfProd(C_auto(s), auto(a))
		CHECK_CALL(glwegadget_half_prod(module, result, automorphism_ksk->enc_s[0], auto_tmp),
		           "half product in automorphism failed");

		// result = -result = -hafProd(C_auto(s), auto(a))
		negate_glwe(module, result, result);

		// auto_tmp = auto_p(b)
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, biv_l, nn, glwe_extract_start_poly(glwe, 1),
		                          glwe_params_l_b(glwe->params), 2 * nn);

		// result += auto_tmp ==> result = -halfProc(c_auto(s), auto(a)) + (0, auto(b))
		pvda_vec_znx_add(module, glwe_extract_start_poly(result, 1), l_b_result, 2 * nn,
		                 glwe_extract_start_poly(result, 1), l_b_result, 2 * nn, auto_tmp, biv_l, nn);

		status = 0;
	cleanup:
		free(auto_tmp);
	}
	else
	{
		PolyBiv* auto_tmp        = new_biv_poly_custom_l(result->params, biv_l);
		GLWECiphertext* glwe_tmp = new_glwe(result->params);

		// auto_tmp = auto_p(a_0)
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, biv_l, nn, glwe_extract_start_poly(glwe, 0),
		                          glwe_params_l_a(glwe->params), (k + 1) * nn);

		// result = halfProd(C_auto(s), auto(a_0))
		CHECK_CALL(glwegadget_half_prod(module, result, automorphism_ksk->enc_s[0], auto_tmp),
		           "half product in automorphism failed");

		for (int i = 1; i < k; ++i)
		{
			// auto_tmp = auto_p(a_i)
			pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, biv_l, nn, glwe_extract_start_poly(glwe, i),
			                          glwe_params_l_a(glwe->params), (k + 1) * nn);

			// result = halfProd(C_auto(s_i), auto(a_i))
			CHECK_CALL_LABEL(glwegadget_half_prod(module, glwe_tmp, automorphism_ksk->enc_s[i], auto_tmp),
			                 "half product in automorphism failed", cleanup2);

			add_glwe(module, result, result, glwe_tmp);
		}

		negate_glwe(module, result, result);
		// auto_tmp = auto_p(b)
		pvda_vec_znx_automorphism(module, automorphism_p, auto_tmp, biv_l, nn, glwe_extract_start_poly(glwe, k),
		                          glwe_params_l_b(glwe->params), (k + 1) * nn);

		// result += auto_tmp ==> result = -sum_i(halfProc(c_auto(s), auto(a_i))) + (0, auto(b))
		pvda_vec_znx_add(module, glwe_extract_start_poly(result, k), l_b_result, (k + 1) * nn,
		                 glwe_extract_start_poly(result, k), l_b_result, (k + 1) * nn, auto_tmp, biv_l, nn);
		status = 0;
	cleanup2:
		free(auto_tmp);
		delete_glwe(glwe_tmp);
	}

	return status;
}
