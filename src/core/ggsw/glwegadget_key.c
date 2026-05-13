#include "glwegadget_key.h"

#include <stdlib.h>

#include "ggsw_params.h"
#include "glwegadget_ciphertext.h"
#include "utils.h"

static inline int allocate_prepared_gadget_array(GLWEGadgetParams* params, GLWEGadgetCiphertextPrep** array, int n)
{
	int i = 0;
	for (int i = 0; i < n; ++i)
	{
		array[i] = new_glwegadget_prep(params);
		CHECK_ALLOC(array[i], "KSK allocation failed");
	}
	return 0;
cleanup:
	for (int i = 0; i < n; ++i)
	{
		delete_glwegadget_prep(array[i]);
	}
	return -1;
}

GLWEAutomorphismKSK* new_automorphism_ksk(GLWEGadgetParams* params)
{
	GLWEAutomorphismKSK* ksk = malloc(sizeof(GLWEAutomorphismKSK));
	CHECK_ALLOC(ksk, "Automorphism KSK allocation failed");
	ksk->params = params;

	ksk->automorphism_p = -2;

	ksk->enc_s = (GLWEGadgetCiphertextPrep**)calloc(params->params_glwe->k, sizeof(GLWEGadgetCiphertextPrep*));
	CHECK_ALLOC(ksk, "Automorphism KSK allocation failed");
	CHECK_CALL(allocate_prepared_gadget_array(params, ksk->enc_s, params->params_glwe->k),
	           "Failed to allocate automorphism KSK key");

	return ksk;
cleanup:
	if (ksk) free((void*)ksk->enc_s);
	free(ksk);
	return NULL;
}

void delete_automorphism_ksk(GLWEAutomorphismKSK* automorphism_ksk)
{
	if (!automorphism_ksk) return;
	for (int i = 0; i < automorphism_ksk->params->params_glwe->k; ++i)
	{
		delete_glwegadget_prep(automorphism_ksk->enc_s[i]);
	}
	free((void*)automorphism_ksk->enc_s);
	free(automorphism_ksk);
}
