#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "spqlios_alias.h"

/**
 * @brief Encrypts the message m into the bivGGSW ciphertext res according to the
 * parameters in the result object.
 *
 * @param module      Additionnal information for backend.
 * @param result      The encrypted message.
 * @param sk_dft      The secret key.
 * @param m_univ      The univariate message (in coefficient space).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_secret_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWESecretKeyDFT* sk_dft,
                        const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a bivGGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate GLWE result ciphertext.
 * @param glwe   The bivariate GLWE input ciphertext.
 * @param ggsw   The bivariate GGSW input ciphertext.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                          const GGSWCiphertext* ggsw);

/**
 * @brief Encrypts message m_univ into GGSW ciphertext res with parameters params
 *
 * TODO: implement
 *
 * @param module The backend module
 * @param result The resulting GGSWCiphertext
 * @param pk     The public key used for encryption
 * @param m_univ The univariate message to encrypt (in coefficient space)
 *
 * */
int ggsw_public_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWEPublicKey* pk, const PolyUniv* m_univ);

#endif  // bivGGSW_H
