#ifndef GLWEGADGET_KEY_H
#define GLWEGADGET_KEY_H

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

#endif
