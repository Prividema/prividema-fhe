#ifndef GGSW_KEY_H
#define GGSW_KEY_H

/**
 * @file ggsw_key.h
 *
 * In this header file, we define the structure representing GGSW secret key in both DFT and iDFT forms.
 */

#include "ggsw_ciphertext.h"

typedef struct ggsw_public_key
{
	uint32_t size;       // number of ciphertexts in public key
	GGSWCiphertext* pk;  // Public key is multiple encryptions of 0
} GGSWPublicKey;

#endif  // GGSW_KEY_H
