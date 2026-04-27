#include "ggsw_ciphertext.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "utils.h"

// bivGGSW Part (begin)

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

// bivGGSW DFT PART (begin)

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
