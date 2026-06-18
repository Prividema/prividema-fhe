#include "glwegadget_ciphertext.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

GLWEGadgetCiphertext* new_glwegadget(const GLWEGadgetParams* params)
{
	GLWEGadgetCiphertext* glwegadget = calloc(sizeof(GLWEGadgetCiphertext), 1);
	CHECK_ALLOC(glwegadget, "alloc failed in GLWEGadget creation");

	glwegadget->params = params;

	glwegadget->mat = calloc(glwegadget_coef_number(params), sizeof(MatBiv));
	CHECK_ALLOC(glwegadget->mat, "alloc failed in GLWEGadget creation");

	return glwegadget;
cleanup:
	free(glwegadget);
	return NULL;
}

void delete_glwegadget(GLWEGadgetCiphertext* glwegadget_ct)
{
	if (!glwegadget_ct) return;
	free(glwegadget_ct->mat);
	free(glwegadget_ct);
}

int print_coefs_gad(const MODULE* module, const GLWEGadgetCiphertext* glwe_gad, const GLWESecretKeyPrepared* sk_prep,
                    int n)
{
	for (int l = 0; l < glwe_gad->params->l_tilde; ++l)
	{
		GLWECiphertext limb = {glwe_gad->params->params_glwe, glwegadget_extract_bivglwe(glwe_gad, l + 1)};
		printf("Gadget lvl %03d: ", l);
		CHECK_CALL(print_coefs_glwe(module, &limb, sk_prep, n, 0), "GLWE printig failed in GLWEGadget printing");
		printf("\n");
	}

	return 0;
cleanup:
	return -1;
}

GLWEGadgetCiphertextPrep* new_glwegadget_prep(const GLWEGadgetParams* params)
{
	GLWEGadgetCiphertextPrep* glwegad_prep = malloc(sizeof(GLWEGadgetCiphertextPrep));
	CHECK_ALLOC(glwegad_prep, "alloc failed in GLWEGadgetPrepared creation");

	glwegad_prep->params = params;

	glwegad_prep->mat = aligned_alloc(64, glwegadget_coef_number(params) * sizeof(MatBivDFT));
	CHECK_ALLOC(glwegad_prep->mat, "alloc failed in GLWEGadgetPrepared creation");
	memset(glwegad_prep->mat, 0, glwegadget_coef_number(params) * sizeof(MatBivDFT));

	return glwegad_prep;
cleanup:
	free(glwegad_prep);
	return NULL;
}

void delete_glwegadget_prep(GLWEGadgetCiphertextPrep* glwegadget_prep_ct)
{
	if (!glwegadget_prep_ct) return;
	free(glwegadget_prep_ct->mat);
	free(glwegadget_prep_ct);
}

VecBiv* glwegadget_extract_bivglwe(GLWEGadgetCiphertext* glwegadget_ct, uint64_t i)
{
	return glwegadget_ct->mat + (i - 1) * glwe_coef_number(glwegadget_ct->params->params_glwe);
}

int glwegadget_secret_encrypt(const MODULE* module, GLWEGadgetCiphertext* result, const GLWESecretKeyPrepared* sk_prep,
                              const PolyUniv* m_univ)
{
	int status = -1;

	const GLWEGadgetParams* params_glwegadget = result->params;
	const GLWEParams* params_glwe             = params_glwegadget->params_glwe;

	uint64_t nn = params_glwe->nn;

	PolyBiv* glwe_biv_msg = new_biv_poly(params_glwe);

	CHECK_ALLOC(glwe_biv_msg, "alloc failed in GLWEGadget encryption");

	for (uint64_t i = 1; i <= params_glwegadget->l_tilde; i++)
	{
		// Computes m_univ / 2^{kappa_tilde*i}
		CHECK_CALL(univ_znx_to_biv(params_glwe, glwe_biv_msg, m_univ, params_glwegadget->kappa_tilde * i),
		           "univ_to_biv failed in compute_phase_ij");

		CHECK_CALL(add_biv_noise(module, params_glwe, glwe_biv_msg, glwe_biv_msg),
		           "Noise addition failed in GLWEGadget encryption");

		// Get the pointer for the result position
		VecBiv* glwe_vec       = glwegadget_extract_bivglwe(result, i);
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec};

		//Compute: bivGLWE(glwe_biv_msg) into glwe_vec (ie, add A * S)
		CHECK_CALL(glwe_secret_encrypt_phase(module, &glwe_ct, sk_prep, glwe_biv_msg),
		           "glwe masking failed in a GLWEGadget encryption");
	}

	status = 0;

cleanup:
	delete_biv(glwe_biv_msg);

	return status;
}

int glwegadget_prepare(const MODULE* module, GLWEGadgetCiphertextPrep* glwegadget_prep_ct,
                       const GLWEGadgetCiphertext* glwegad_ct)
{
	int status = -1;

	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	CHECK_CALL(pvda_vmp_prepare_contiguous(module, glwegadget_prep_ct->mat, glwegad_ct->mat, nrows, ncols),
	           "VMP prepare for GLWEGadget prepare failed");

	status = 0;
cleanup:
	return status;
}

int glwegadget_packed_secret_encrypt(const MODULE* module, GLWECiphertext* result,
                                     const GLWEGadgetParams* params_glwegad, const GLWESecretKeyPrepared* sk_prep,
                                     const PolyUniv* m_univ, uint64_t d)
{
	int status                    = -1;
	const GLWEParams* params_glwe = result->params;

	uint64_t nn = params_glwe->nn;

	PolyBiv* glwe_biv_msg = new_biv_poly(params_glwe);
	PolyBiv* glwe_biv_tmp = new_biv_poly(params_glwe);
	PolyUniv* pol_encrypt = new_univ(params_glwe);
	CHECK_ALLOC(glwe_biv_msg, "alloc failed in GLWEGadget encryption");
	CHECK_ALLOC(glwe_biv_tmp, "alloc failed in GLWEGadget encryption");
	CHECK_ALLOC(pol_encrypt, "alloc failed in GLWEGadget packed encryption");

	assert(params_glwegad->l_tilde * d <= nn);

	int64_t divlog = next_pow2_log(d * params_glwegad->l_tilde);
	for (uint64_t i = 1; i <= params_glwegad->l_tilde; i++)
	{
		//Can be optimised
		memset(pol_encrypt, 0, poly_univ_bytes(params_glwe));
		memcpy(pol_encrypt + (i - 1) * d, m_univ, d * sizeof(PolyUniv));
		CHECK_CALL(univ_znx_to_biv(params_glwe, glwe_biv_tmp, pol_encrypt, params_glwegad->kappa_tilde * i + divlog),
		           "univ_to_biv failed in compute_phase_ij");
		add_biv_poly(module, params_glwe, glwe_biv_msg, glwe_biv_msg, glwe_biv_tmp);
	}
	CHECK_CALL(add_biv_noise(module, params_glwe, glwe_biv_msg, glwe_biv_msg),
	           "Noise addition failed in GLWEGadget encryption");

	CHECK_CALL(glwe_secret_encrypt_phase(module, result, sk_prep, glwe_biv_msg),
	           "glwe masking failed in a GLWEGadget encryption");

	status = 0;

cleanup:
	delete_biv(glwe_biv_msg);
	delete_biv(glwe_biv_tmp);
	free(pol_encrypt);

	return status;
}
