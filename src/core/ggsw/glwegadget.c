#include "glwegadget.h"

#include <math.h>
#include <stdint.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "glwegadget_ciphertext.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

VecBiv* glwegadget_extract_bivglwe(GLWEGadgetCiphertext* glwegadget_ct, uint64_t i)
{
	// bivGLWE parameters
	const GLWEParams* params_glwe = glwegadget_ct->params->params_glwe;

	return glwegadget_ct->mat + (i - 1) * glwe_coef_number(params_glwe);
}

int glwegadget_secret_encrypt(const MODULE* module, GLWEGadgetCiphertext* result, const GLWESecretKeyDFT* sk_dft,
                              const PolyUniv* m_univ)
{
	const GLWEGadgetParams* params_glwegadget = result->params;

	int status = -1;

	const GLWEParams* params_glwe = params_glwegadget->params_glwe;

	// bivGLWE parameters
	uint64_t nn = params_glwe->nn;
	uint64_t k  = params_glwe->k;

	PolyUnivRnX* tmp_sp1  = new_univ_rnx(params_glwe);
	PolyBiv* glwe_biv_msg = new_biv_poly(params_glwe);

	CHECK_ALLOC(tmp_sp1, "alloc failed in GLWEGadget encryption");
	CHECK_ALLOC(glwe_biv_msg, "alloc failed in GLWEGadget encryption");

	for (uint64_t i = 1; i <= params_glwegadget->l_tilde; i++)
	{
		// Computes m_univ / 2^{kappa_tilde*i}
		for (uint64_t p = 0; p < nn; p++) tmp_sp1[p] = ldexp((double)m_univ[p], -params_glwegadget->kappa_tilde * i);

		// Adds the error
		CHECK_CALL(add_normal_random_vec(tmp_sp1, nn, tmp_sp1, 0.0, params_glwe->sigma),
		           "Error addition failed in GLWEGadget encryption");
		// Compute the base-2^kappa decomposition
		CHECK_CALL(univ_rnx_to_biv(params_glwe, glwe_biv_msg, tmp_sp1), "univ_to_biv failed in compute_phase_ij");

		// Get the pointer for the result position
		VecBiv* glwe_vec       = glwegadget_extract_bivglwe(result, i);
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec};

		//Compute: bivGLWE(glwe_biv_msg) into glwe_vec (ie, add A * S)
		CHECK_CALL(glwe_secret_encrypt_phase(module, &glwe_ct, sk_dft, glwe_biv_msg),
		           "glwe masking failed in a GLWEGadget encryption");
	}

	status = 0;

cleanup:
	free(glwe_biv_msg);
	delete_univ_rnx(tmp_sp1);

	return status;
}

int glwegadget_prepare(const MODULE* module, GLWEGadgetCiphertextPrep* glwegadget_prep_ct,
                       const GLWEGadgetCiphertext* glwegad_ct)
{
	int status = -1;

	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	uint64_t nn  = glwegadget_prep_ct->params->params_glwe->nn;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	CHECK_CALL(pvda_vmp_prepare_contiguous(module, glwegadget_prep_ct->mat, glwegad_ct->mat, nrows, ncols),
	           "VMP prepare for GLWEGadget prepare failed");

	status = 0;
cleanup:
	return status;
}

int glwegadget_half_prod(const MODULE* module, GLWECiphertext* result,
                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBiv* a)
{
	int status = -1;

	//TODO: assert size compatibility
	size_t nrows = glwegadget_prep_ct->params->l_tilde;
	uint64_t nn  = glwegadget_prep_ct->params->params_glwe->nn;
	size_t ncols = glwe_params_n_limbs(glwegadget_prep_ct->params->params_glwe);

	GLWECiphertextDFT* glwe_dft = new_glwe_dft(result->params);

	CHECK_CALL(pvda_vmp_apply_dft(module, glwe_dft->vec, ncols, a, nrows, nn, glwegadget_prep_ct->mat, nrows, ncols),
	           "vmp apply falied in half product");

	CHECK_CALL(glwe_dft_to_coef(module, result, glwe_dft),
	           "conversion from GLWE DFT to coefs failed in GLWEGadget half product");

cleanup:
	delete_glwe_dft(glwe_dft);

	return status;
}
