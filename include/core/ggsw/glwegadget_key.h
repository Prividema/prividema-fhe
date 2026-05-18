#ifndef GLWEGADGET_KEY_H
#define GLWEGADGET_KEY_H

#include <stdint.h>

#include "glwegadget_ciphertext.h"

/**
 * @file glwegadget_key.h
 *
 * @brief Contains data structures related GLWEGadget specific keys, mainly key-switching-keys.
 *
 * So far only KSK for automprphisms are here
 */

typedef struct glwe_automorphism_keyswitching_key_t
{
	int automorphism_p;  // -2 represents uninitialized value
	GLWEGadgetParams* params;
	GLWEGadgetCiphertextPrep** enc_s;
} GLWEAutomorphismKSK;

GLWEAutomorphismKSK* new_automorphism_ksk(GLWEGadgetParams* params);

void delete_automorphism_ksk(GLWEAutomorphismKSK* automorphism_ksk);

typedef struct glwe_automorphism_keyswitching_key_collection_t
{
	uint64_t size;
	GLWEAutomorphismKSK** keys;
} GLWEAutomorphismKSKCollection;

GLWEAutomorphismKSKCollection* new_automorphism_ksk_collection(uint64_t size);

GLWEAutomorphismKSK* glwegadget_ksk_collection_put_key(GLWEAutomorphismKSKCollection* collection,
                                                       GLWEAutomorphismKSK* key, uint64_t pos);

GLWEAutomorphismKSK* glwegadget_ksk_collection_get_key(const GLWEAutomorphismKSKCollection* collection, uint64_t pos);

uint64_t delete_automorphism_ksk_collection(GLWEAutomorphismKSKCollection* automorphism_ksk, int deallocate_ksks);

#endif
