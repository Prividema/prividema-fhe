#include "glwegadget_key.h"

#include <stdint.h>
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

GLWEAutomorphismKSKCollection* new_automorphism_ksk_collection(uint64_t size)
{
	GLWEAutomorphismKSKCollection* collection = malloc(sizeof(GLWEAutomorphismKSKCollection));
	CHECK_ALLOC(collection, "Failed allocation of KSK collection");
	collection->keys = calloc(size, sizeof(GLWEAutomorphismKSK*));
	CHECK_ALLOC(collection->keys, "Failed allocation of KSK collection");
	collection->size = size;
	return collection;

cleanup:

	if (collection) free(collection->keys);
	free(collection);
	return NULL;
}

GLWEAutomorphismKSK* glwegadget_ksk_collection_put_key(GLWEAutomorphismKSKCollection* collection,
                                                       GLWEAutomorphismKSK* key, uint64_t pos)
{
	if (pos >= collection->size)
	{
		return (void*)-1;
	}
	GLWEAutomorphismKSK* tmp = collection->keys[pos];
	collection->keys[pos]    = key;
	return tmp;
}

GLWEAutomorphismKSK* glwegadget_ksk_collection_get_key(const GLWEAutomorphismKSKCollection* collection, uint64_t pos)
{
	if (pos >= collection->size) return NULL;
	return collection->keys[pos];
}

uint64_t delete_automorphism_ksk_collection(GLWEAutomorphismKSKCollection* automorphism_ksk, int deallocate_ksks)
{
	uint64_t count = 0;
	if (deallocate_ksks)
	{
		for (uint64_t i = 0; i < automorphism_ksk->size; ++i)
		{
			if (automorphism_ksk->keys[i]) ++count;
			delete_automorphism_ksk(automorphism_ksk->keys[i]);
		}
	}
	free(automorphism_ksk->keys);
	free(automorphism_ksk);
	return count;
}
