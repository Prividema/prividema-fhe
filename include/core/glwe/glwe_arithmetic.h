#ifndef bivGLWE_H
#define bivGLWE_H

#include "glwe_ciphertext.h"
#include "glwe_key.h"

/**
 * @brief Normalizes a bivGLWE ciphertext.
 *
 * @param module Additionnal information for backend.
 * @param res The result normalized bivGLWE ciphertext.
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int normalize_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* glwe);

/**
 * @brief Adds two bivariate GLWE ciphertexts.
 *

 * @param module Additionnal information for backend.
 * @param res The result GLWE ciphertext.
 * @param glwe_lhs The left-hand side GLWE ciphertext.
 * @param glwe_rhs The right-hand side GLWE ciphertext.
 */
void add_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs);

/**
 * @brief Subtracts two bivariate GLWE ciphertexts.
 *
 * @param module Additionnal information for backend.
 * @param res The result GLWE ciphertext.
 * @param glwe_lhs The left-hand side GLWE ciphertext.
 * @param glwe_rhs The right-hand side GLWE ciphertext.
 */
void sub_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* glwe_lhs,
              const GLWECiphertext* glwe_rhs);

/**
 * @brief Negates (inverts the sign of) a GLWE
 *
 * @param module Additionnal information for backend.
 * @param res The result GLWE ciphertext (can be the same as the input for in-place negation).
 * @param glwe The GLWE ciphertext to invert.
 */
void negate_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* glwe);

/**
 * @brief Multiply a bivGLWE ciphertext by a \ZnX polynomial.
 *
 * @param module Additionnal information for backend.
 * @param res The result GLWE ciphertext.
 * @param u The \ZnX polynomial.
 * @param glwe The GLWE ciphertext.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int const_mult_glwe(const MODULE* module, GLWECiphertext* res, const PolyUnivDFT* u, const GLWECiphertext* glwe);

/**
 * @brief Adds two bivGLWE ciphertexts.
 *
 * @param res_dft The result bivGLWE ciphertext in the DFT domain.
 * @param glwe_lhs_dft The left-hand side bivGLWE ciphertext in the DFT domain.
 * @param glwe_rhs_dft The right-hand side bivGLWE ciphertext in the DFT domain.
 */
void add_glwe_dft(GLWECiphertextDFT* res_dft, const GLWECiphertextDFT* glwe_lhs_dft,
                  const GLWECiphertextDFT* glwe_rhs_dft);

/**
 * @brief Multiply a bivGLWE ciphertext by a \ZnX polynomial in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param res_dft The result bivGLWE ciphertext in the DFT domain.
 * @param u The \ZnX polynomial.
 * @param glwe_dft The bivGLWE ciphertext in the DFT domain.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int const_mult_glwe_dft(const MODULE* module, GLWECiphertextDFT* res_dft, const PolyUnivDFT* u,
                        const GLWECiphertextDFT* glwe_dft);

/**
 * @brief Compute the polynomial product of c and d, component-wise in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param res_dft The result in the DFT domain.
 * @param res_size The result's size.
 * @param c_dft The left-hand side polynomial in the DFT domain .
 * @param c_size The left-hand size of c_dft.
 * @param d_dft The right-hand side polynomial in the DFT domain.
 * @param d_size The right-hand size of c_dft.
 *
 * @remark `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the
 * size of res_dft.
 */
void mult_vec_znx_dft(const MODULE* module, double* res_dft, int64_t res_size, const double* c_dft, int64_t c_size,
                      const double* d_dft, int64_t d_size);

#endif
