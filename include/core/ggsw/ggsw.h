#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_ciphertext.h"
#include "glwe_key.h"
#include "spqlios_alias.h"

/**
 * @brief Encrypts the message m into bivGGSW ciphertext res with parameters params.
 *
 * @param module      Additionnal information for backend.
 * @param params_ggsw The encryption params.
 * @param result      The encrypted message.
 * @param sk_dft      The secret key.
 * @param m_univ      The message.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise.
 */
int ggsw_secret_encrypt(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertext* result,
                        const GLWESecretKeyDFT* sk_dft, const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a bivGGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate GLWE result ciphertext.
 * @param glwe   The bivariate GLWE input ciphertext.
 * @param ggsw   The bivariate GGSW input ciphertext.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise.
 */
int ggsw_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                          const GGSWCiphertext* ggsw);

// The three functions below are not implemented

// GGSWPublicKey is a struct encapsulating everything regarding the public
// key.
/* Encrypts message m into GGSW ciphertext res with parameters params */
void ggsw_public_encrypt(GGSWCiphertext* result,     // result
                         const GLWEPublicKey* pk,    // secret key
                         const int64_t* m,           // message
                         const GGSWCtParams* params  // parameters
);

// PartialGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_secret_encrypt(PartialGGSWCiphertext* result,  // result
                             const GLWESecretKey* sk,        // secret key
                             const int64_t* m,               // message
                             const GGSWCtParams* params      // parameters
);

/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_public_encrypt(PartialGGSWCiphertext* result,  // result
                             const GLWEPublicKey* pk,        // public key
                             const int64_t* m,               // message
                             const GGSWCtParams* params      // parameters
);

//! GGSW IN DFT PART (begin)

/**
 * @brief Encrypts the message m into bivGGSW ciphertext res with parameters params in the DFT domain.
 *
 * @param module      Additionnal information for backend.
 * @param params_ggsw The encryption params.
 * @param res_dft     The encrypted message.
 * @param sk_dft      The secret key.
 * @param m_univ      The message.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise.
 */
int ggsw_secret_encrypt_dft(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertextDFT* result_dft,
                            const GLWESecretKeyDFT* sk_dft, const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a biv bivGGSW.
 *
 * @param module   Additionnal information for backend.
 * @param res_dft  The bivariate GLWE result ciphertext in the DFT domain.
 * @param glwe_dft The bivariate GLWE input ciphertext in the DFT domain.
 * @param ggsw_dft The bivariate GGSW input ciphertext in the DFT domain.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise.
 */
int ggsw_external_product_dft(const MODULE* module, GLWECiphertextDFT* result_dft, const GLWECiphertextDFT* glwe_dft,
                              const GGSWCiphertextDFT* ggsw_dft);

#endif  // bivGGSW_H
