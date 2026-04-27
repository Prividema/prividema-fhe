#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "spqlios_alias.h"

/**
 * @brief Normalizes a GGSW ciphertext.
 *
 * @param module Additionnal information for backend.
 * @param result The oupput GGSWCiphertext. It can be the same as the output for in-place normalization
 * @param ggsw   The input GGSWCiphertext
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw);

/**
 * @brief Adds two GGSW ciphertexts that have the same parameters.
 *
 * @param result    The resulting GGSWCiphertext
 * @param ggsw_lhs  One of the input GGSWCiphertext
 * @param ggsw_rhs  The other input GGSWCiphertext
 */
void add_ggsw(GGSWCiphertext* result, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs);

/**
 * @brief Multiplies a GGSW ciphertext by a constant in \ZnX
 *
 * @param module    Additionnal information for backend.
 * @param result    The resulting GGSWCiphertext
 * @param ggsw      The input GGSWCiphertext.
 * @param cte       The constant in \ZnX
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* cte);

/**
 * @brief Adds two GGSW ciphertexts in the DFT space  that have the same parameters
 *
 * @param result_dft   The resulting GGSWCiphertextDFT
 * @param ggsw_lhs_dft One of the input GGSWCiphertextDFT.
 * @param ggsw_rhs_dft The other GGSWCiphertextDFT.
 */
void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft,
                  const GGSWCiphertextDFT* ggsw_rhs_dft);

/**
 * @brief Multiplies a GGSW ciphertext in the DFT space by a polynomial from \ZnX also in DFT space
 *
 * @param module    Additionnal information for backend.
 * @param result_dft The output GGSWCiphertextDFT
 * @param ggsw_dft      The input GGSWCiphertextDFT.
 * @param u_dft       The \ZnX polynomial in the DFT space.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft);

/**
 * @brief Computes the external product between a GLWE and a GGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The GLWE result ciphertext.
 * @param glwe   The GLWE input ciphertext.
 * @param ggsw   The GGSW input ciphertext.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                          const GGSWCiphertext* ggsw);

#endif  // bivGGSW_H
