#include "glwegad_ciphertext.h"

#include "ggsw_params.h"
#include "glwe_ciphertext.h"
#include "rng.h"
#include "utils.h"

size_t glwegad_coef_number(const GLWEGadParams* params_glwegad)
{
	return params_glwegad->l_tilde * glwe_coef_number(params_glwegad->params_glwe);
}

GLWEGadCiphertext* new_glwegad(const GLWEGadParams* params)
{
	GLWEGadCiphertext* glwegad = malloc(sizeof(GLWEGadCiphertext));
	CHECK_ALLOC(glwegad, "malloc failed in partial GGSW creation");

	glwegad->params = params;

	glwegad->mat = calloc(glwegad_coef_number(params), sizeof(MatBiv));
	CHECK_ALLOC(glwegad, "malloc failed in partial GGSW creation");

	return glwegad;
cleanup:
	free(glwegad);
	return NULL;
}

void delete_glwegad(GLWEGadCiphertext* glwegad_ct)
{
	if (!glwegad_ct) return;
	free(glwegad_ct->mat);
	free(glwegad_ct);
}
