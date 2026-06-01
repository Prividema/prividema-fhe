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
	int automorphism_p;                ///< -2 represents uninitialized value, 0 if the KSK is not for an automorphism
	GLWEGadgetParams* params;          ///< The GLWEGadgetParams associated with the gadgets used for key switching
	GLWEGadgetCiphertextPrep** enc_s;  ///< A vector with pointers to the k prepared ciphertetxts of f(sk_i)
} GLWEAutomorphismKSK;

/**
 * @brief Creates a new Key-switching-key (KSK)
 *
 * @param params The GLWEGadgetParams associated with the KSK
 *
 * @return NULL in case of failure, a new KSK otherwise
 *
 */
GLWEAutomorphismKSK* new_automorphism_ksk(GLWEGadgetParams* params);

/**
 * @brief Deletes a KSK
 *
 * @param automorphism_ksk The KSK to delete
 */
void delete_automorphism_ksk(GLWEAutomorphismKSK* automorphism_ksk);

#endif
