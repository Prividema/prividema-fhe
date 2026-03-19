#ifndef PARTIAL_GGSW_H
#define PARTIAL_GGSW_H

#include "bivariate_polynomial.h"
#include "glwe_key.h"
#include "glwegad_ciphertext.h"

/**
 * @brief Encrypts a univariate integer pollynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 * @params module The underlying compute module
 * @params result The resulting GLWEGadget ciphertext
 * @params sk_dft A prepared secret key
 * @params m_univ The univariate plaintext
 *
 */
int glwegad_secret_encrypt(const MODULE* module,
                           GLWEGadCiphertext* result,       // result
                           const GLWESecretKeyDFT* sk_dft,  // secret key
                           const PolyUniv* m_univ           // message
);

/**
 * @brief Encrypts a univariate integer pollynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 *  TODO: implement
 *
 * @params module The underlying compute module
 * @params result The resulting GLWEGadget ciphertext
 * @params pk     A public key
 * @params m_univ The univariate plaintext
 *
 */
void glwegad_public_encrypt(const MODULE* module,
                            GLWEGadCiphertext* result,  // result
                            const GLWEPublicKey* pk,    // public key
                            const PolyUniv* m_univ      // message
);

/**
 * @brief Retreive a bivglwe from a GLWEGadget (halfGGSW for k = 1)
 *
 * @param glwegad_ct     The GLWEGadget ciphertext
 * @param i            The index of the GLWE to retrievie, from 1 to l_tilde included
 *
 */
VecBiv* glwegad_retrieve_bivglwe(GLWEGadCiphertext* glwegad_ct, int64_t i);

#endif  // !DEBUG
