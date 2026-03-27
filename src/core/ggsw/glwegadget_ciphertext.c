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
