#ifndef GGSW_KEY_H
#define GGSW_KEY_H

/**
 * @file ggsw_key.h
 *
 * In this header file, we define the structure representing GGSW public keys.
 * Secret keys are the same for GLWE and GGSW
 */

#include "ggsw_ciphertext.h"

/**
 * @brief A GGSW public key
 *
 * In other words, a collection of encryptions of 0 under the same secret key.
 *
 */
typedef struct ggsw_public_key
{
	uint32_t size;       ///< number of ciphertexts in public key
	GGSWCiphertext* pk;  ///< Public key is multiple encryptions of 0
} GGSWPublicKey;

#endif  // GGSW_KEY_H
