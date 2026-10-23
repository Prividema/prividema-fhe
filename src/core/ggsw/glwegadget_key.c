#include "glwegadget_key.h"

#include <stdint.h>
#include <stdlib.h>

#include "ggsw_params.h"
#include "glwegadget_ciphertext.h"
#include "utils.h"

static inline int allocate_prepared_gadget_array(const PvdaBackend* module, const GLWEGadgetParams* params_glwegadget,
                                                 GLWEGadgetCiphertextPrep** array, int n)
{
	int i = 0;
	for (int i = 0; i < n; ++i)
	{
		array[i] = new_glwegadget_prep(module, params_glwegadget);
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

GLWEKSK* new_glwe_ksk(const PvdaBackend* module, const GLWEGadgetParams* params_glwegadget)
{
	GLWEKSK* glwe_ksk = malloc(sizeof(GLWEKSK));
	CHECK_ALLOC(glwe_ksk, "KSK allocation failed");
	glwe_ksk->params = params_glwegadget;

	glwe_ksk->enc_s =
	    (GLWEGadgetCiphertextPrep**)calloc(params_glwegadget->params_glwe->k, sizeof(GLWEGadgetCiphertextPrep*));
	CHECK_ALLOC(glwe_ksk->enc_s, "KSK gadget array allocation failed");
	CHECK_CALL(
	    allocate_prepared_gadget_array(module, params_glwegadget, glwe_ksk->enc_s, params_glwegadget->params_glwe->k),
	    "Failed to allocate the KSK prepared gadgets");

	return glwe_ksk;
cleanup:
	if (glwe_ksk) free((void*)glwe_ksk->enc_s);
	free(glwe_ksk);
	return NULL;
}

void delete_glwe_ksk(GLWEKSK* glwe_ksk)
{
	if (!glwe_ksk) return;
	for (int i = 0; i < glwe_ksk->params->params_glwe->k; ++i)
	{
		delete_glwegadget_prep(glwe_ksk->enc_s[i]);
	}
	free(glwe_ksk->enc_s);
	free(glwe_ksk);
}

GLWEAutomorphismKey* new_automorphism_key(const PvdaBackend* module, const GLWEGadgetParams* params_glwegadget)
{
	GLWEAutomorphismKey* automorphism_key = malloc(sizeof(GLWEAutomorphismKey));
	CHECK_ALLOC(automorphism_key, "Automorphism key allocation failed");

	automorphism_key->automorphism_p = -2;

	automorphism_key->glwe_ksk = new_glwe_ksk(module, params_glwegadget);
	CHECK_ALLOC(automorphism_key->glwe_ksk, "Automorphism key allocation failed");

	return automorphism_key;
cleanup:
	if (automorphism_key) delete_glwe_ksk(automorphism_key->glwe_ksk);
	free(automorphism_key);
	return NULL;
}

void delete_automorphism_key(GLWEAutomorphismKey* automorphism_key)
{
	if (!automorphism_key) return;
	delete_glwe_ksk(automorphism_key->glwe_ksk);
	free(automorphism_key);
}

GLWEAutomorphismKeyCollection* new_automorphism_key_collection(uint64_t size)
{
	GLWEAutomorphismKeyCollection* collection = malloc(sizeof(GLWEAutomorphismKeyCollection));
	CHECK_ALLOC(collection, "Failed allocation of automorphism key collection");
	collection->keys = calloc(size, sizeof(GLWEAutomorphismKey*));
	CHECK_ALLOC(collection->keys, "Failed allocation of automorphism key collection");
	collection->size = size;
	return collection;

cleanup:

	if (collection) free(collection->keys);
	free(collection);
	return NULL;
}

GLWEAutomorphismKey* glwegadget_key_collection_put_key(GLWEAutomorphismKeyCollection* collection,
                                                       GLWEAutomorphismKey* automorphism_key, uint64_t pos)
{
	if (pos >= collection->size)
	{
		return (void*)-1;
	}
	GLWEAutomorphismKey* tmp = collection->keys[pos];
	collection->keys[pos]    = automorphism_key;
	return tmp;
}

GLWEAutomorphismKey* glwegadget_key_collection_get_key(const GLWEAutomorphismKeyCollection* collection, uint64_t pos)
{
	if (pos >= collection->size) return NULL;
	return collection->keys[pos];
}

uint64_t delete_automorphism_key_collection(GLWEAutomorphismKeyCollection* automorphism_keys, int deallocate_keys)
{
	uint64_t count = 0;
	if (deallocate_keys)
	{
		for (uint64_t i = 0; i < automorphism_keys->size; ++i)
		{
			if (automorphism_keys->keys[i]) ++count;
			delete_automorphism_key(automorphism_keys->keys[i]);
		}
	}
	free(automorphism_keys->keys);
	free(automorphism_keys);
	return count;
}
