#ifndef PARTIAL_GGSW_H
#define PARTIAL_GGSW_H

#include "bivariate_polynomial.h"
#include "glwe_key.h"
#include "glwegadget_ciphertext.h"

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
int glwegadget_secret_encrypt(const MODULE* module,
                              GLWEGadgetCiphertext* result,    // result
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
void glwegadget_public_encrypt(const MODULE* module,
                               GLWEGadgetCiphertext* result,  // result
                               const GLWEPublicKey* pk,       // public key
                               const PolyUniv* m_univ         // message
);

/**
 * @brief Retreive a bivglwe from a GLWEGadget (halfGGSW for k = 1)
 *
 * @param glwegadget_ct     The GLWEGadget ciphertext
 * @param i            The index of the GLWE to retrievie, from 1 to l_tilde included
 *
 */
VecBiv* glwegadget_extract_bivglwe(GLWEGadgetCiphertext* glwegadget_ct, uint64_t i);

/**
 *
 * @brief Computes the half-external product between a GLWEGadget and
 * a GLWE ciphertext
 *
 */
int glwegadget_half_prod(const MODULE* module, GLWEGadgetCiphertext* result,
                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBiv* a);

#endif  // !DEBUG
