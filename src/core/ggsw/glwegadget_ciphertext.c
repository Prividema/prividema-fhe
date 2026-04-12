#include "glwegadget_ciphertext.h"

#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "rng.h"
#include "utils.h"

size_t glwegadget_coef_number(const GLWEGadgetParams* params_glwegadget)
{
	return params_glwegadget->l_tilde * glwe_coef_number(params_glwegadget->params_glwe);
}

GLWEGadgetCiphertext* new_glwegadget(const GLWEGadgetParams* params)
{
	GLWEGadgetCiphertext* glwegadget = malloc(sizeof(GLWEGadgetCiphertext));
	CHECK_ALLOC(glwegadget, "alloc failed in GLWEGadget creation");

	glwegadget->params = params;

	glwegadget->mat = calloc(glwegadget_coef_number(params), sizeof(MatBiv));
	CHECK_ALLOC(glwegadget, "alloc failed in GLWEGadget creation");

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

GLWEGadgetCiphertextPrep* new_glwegadget_prep(const GLWEGadgetParams* params)
{
	GLWEGadgetCiphertextPrep* glwegad_prep = malloc(sizeof(GLWEGadgetCiphertextPrep));
	CHECK_ALLOC(glwegad_prep, "alloc failed in GLWEGadgetPrepared creation");

	glwegad_prep->params = params;

	glwegad_prep->mat = calloc(glwegadget_coef_number(params), sizeof(MatBivDFT));
	CHECK_ALLOC(glwegad_prep->mat, "alloc failed in GLWEGadgetPrepared creation");

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
