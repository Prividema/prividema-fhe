#ifndef GGSW_KEY_H
#define GGSW_KEY_H

/**
 * @file ggsw_key.h
 *
 * In this header file, we define the structure representing GGSW public keys.
 * Secret keys are the same for GLWE and GGSW
 */

#include "ggsw_ciphertext.h"
#include "glwe_key.h"
#include "spqlios_alias.h"

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

/**
 * @brief Generates a KSK for converting GLWEGadgets into GGSWs (a GGSW(-s))
 *
 * @param module the backend module
 * @param ggsw_ksks   The destination set of GGSWs. This should be an array of GGSWCiphertextPrep pointers of size k.
 *                    For values set to NULL, this function will allocate them and transfer responsibility to the caller.
 *                    Otherwise, it will try to use the existing GGSWCiphertextPreps
 * @param ggsw_params GGSW parameters used for the destination and intermetdiate results
 * @param sk_prep     The secret key to use to prepare the KSK
 */
int generate_glwegad_to_ggsw_ksk(const MODULE* module, GGSWCiphertextPrep** ggsw_ksks, const GGSWParams* ggsw_params,
                                 const GLWESecretKeyPrepared* sk_prep);

#endif  // GGSW_KEY_H
