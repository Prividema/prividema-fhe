#ifndef PARTIAL_GGSW_H
#define PARTIAL_GGSW_H

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwegadget_ciphertext.h"

/**
 * @brief Computes the half-external product between a GLWEGadget and
 * a GLWE ciphertext
 *
 * @param module             The backend module
 * @param result             The resulting GLWE ciphertext
 * @param glwegadget_prep_ct A prepared GLWEGadget to multiply with
 * @param a                  An integer bivariate polynomial
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_half_prod(const MODULE* module, GLWECiphertext* result,
                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBiv* a);

/**
 * @brief Computes the half-external product between a GLWEGadget and
 * a GLWE ciphertext, with inputs and outputs in the DFT domain
 *
 * @param module                The backend module
 * @param result_dft            The resulting GLWE ciphertext in DFT domain
 * @param glwegadget_prep_ct    A prepared GLWEGadget to multiply with
 * @param a_dft                 An bivariate polynomial in the DFT domain
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_half_prod_dft_to_dft(const MODULE* module, GLWECiphertextDFT* result_dft,
                                    const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBivDFT* a_dft);

#endif  // !DEBUG
