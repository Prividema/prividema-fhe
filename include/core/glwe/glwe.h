#ifndef bivGLWE_H
#define bivGLWE_H

#include "glwe_key.h"

//! bivGLWE PART (begin)

/**
 * @brief Encrypts a phase (message + noise) and puts it in result.
 *
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase message + noise.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval  -`0` othwerwise.
 */
int glwe_secret_encrypt_phase(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                              const PolyBiv* phase);

/**
 * @brief Encrypts a univariate message
 *
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ_rnx  The univariate secret message to encrypt
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval  -`0` othwerwise.
 */
int glwe_secret_encrypt_rnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                            const PolyUnivRnX* m_univ_rnx);

/**
 * @brief Encrypts a univariate message
 *
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ_tnx  The univariate secret message to encrypt
 *
 * @retval -1 if an error occurs.
 * @retval  -0 othwerwise.
 */
int glwe_secret_encrypt_tnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                            const PolyUnivTnX* m_univ_tnx);

/**
 * @brief Demasks the ciphertext into the
 * bivariate phase (message + noise) and puts it in result.
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate phase.
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval  -`0` othwerwise.
 */
int glwe_secret_decrypt(const MODULE* module, PolyBiv* result, const GLWESecretKeyDFT* sk_dft,
                        const GLWECiphertext* glwe);

#endif
