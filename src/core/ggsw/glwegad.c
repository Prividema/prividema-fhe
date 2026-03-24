#include "glwegad.h"

#include <math.h>
#include <stdint.h>

#include "glwe.h"
#include "glwe_ciphertext.h"
#include "glwegad_ciphertext.h"
#include "rng.h"
#include "utils.h"

VecBiv* glwegad_retrieve_bivglwe(GLWEGadCiphertext* glwegad_ct, int64_t i)
{
	// bivGLWE parameters
	const GLWEParams* params_glwe = glwegad_ct->params->params_glwe;

	return glwegad_ct->mat + (i - 1) * glwe_coef_number(params_glwe);
}

int glwegad_secret_encrypt(const MODULE* module, GLWEGadCiphertext* result, const GLWESecretKeyDFT* sk_dft,
                           const PolyUniv* m_univ)
{
	const GLWEGadParams* params_glwegad = result->params;

	int status = -1;

	const GLWEParams* params_glwe = params_glwegad->params_glwe;

	// bivGLWE parameters
	uint64_t nn = params_glwe->nn;
	uint64_t k  = params_glwe->k;

	// TODO: find out why it is being set to poly_univ_bytes even though it is a different data type
	double* tmp_sp1 = NULL;

	PolyBiv* glwe_biv_msg = NULL;

	tmp_sp1 = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(tmp_sp1, "malloc failed in ggsw_secret_encrypt");
	glwe_biv_msg = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg, "malloc failed in ggsw_secret_encrypt");

	for (uint64_t i = 1; i <= params_glwegad->l_tilde; i++)
	{
		// Computes m_univ / 2^{kappa_tilde*i}
		for (uint64_t p = 0; p < nn; p++) tmp_sp1[p] = ldexp((double)m_univ[p], -params_glwegad->kappa_tilde * i);

		// Adds the error
		CHECK_CALL(add_normal_random_vec(tmp_sp1, nn, tmp_sp1, 0.0, params_glwe->sigma),
		           "Error addition failed in partialGGSW encryption");
		// Compute the base-2^kappa decomposition
		CHECK_CALL(univ_to_biv(params_glwe, glwe_biv_msg, tmp_sp1), "univ_to_biv failed in compute_phase_ij");

		// Get the pointer for the result position
		VecBiv* glwe_vec       = glwegad_retrieve_bivglwe(result, i);
		GLWECiphertext glwe_ct = {params_glwe, glwe_vec};

		//Compute: bivGLWE(glwe_biv_msg) into glwe_vec (ie, add A * S)
		CHECK_CALL(glwe_secret_masking(module, &glwe_ct, sk_dft, glwe_biv_msg),
		           "glwe_secret_masking_ggsw_lib failed in ggsw_secret_encrypt");
	}

	status = 0;

cleanup:
	free(glwe_biv_msg);
	free(tmp_sp1);

	return status;
}

int glwegad_half_prod(const MODULE* module, GLWEGadCiphertext* result, const GLWEGadCiphertextPrep* glwegad_prep_ct,
                      const PolyBiv* a)
{
	int status = -1;

	//TODO: assert size compatibility
	size_t nrows = glwegad_prep_ct->params->l_tilde;
	uint64_t nn  = glwegad_prep_ct->params->params_glwe->nn;
	size_t ncols = glwegad_prep_ct->params->params_glwe->n_limbs;

	GLWECiphertextDFT* glwe_dft = new_glwe_dft(result->params->params_glwe);

	CHECK_CALL(pvda_vmp_apply_dft(module, glwe_dft->vec, ncols, a, nrows, nn, glwegad_prep_ct->mat, nrows, ncols),
	           "vmp apply falied in half product");

	CHECK_CALL(pvda_vec_znx_idft(module, result->mat, ncols, glwe_dft->vec, ncols), "iDFT failed in half product");

cleanup:
	delete_glwe_dft(glwe_dft);

	return status;
}
