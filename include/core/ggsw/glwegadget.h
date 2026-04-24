#ifndef PARTIAL_GGSW_H
#define PARTIAL_GGSW_H

#include "bivariate_polynomial.h"
#include "glwe_key.h"
#include "glwegadget_ciphertext.h"

/**
 * @brief Encrypts a univariate integer polynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 * @param module The underlying compute module
 * @param result The resulting GLWEGadget ciphertext
 * @param sk_dft A prepared secret key
 * @param m_univ The univariate plaintext
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_secret_encrypt(const MODULE* module, GLWEGadgetCiphertext* result, const GLWESecretKeyDFT* sk_dft,
                              const PolyUniv* m_univ);

/**
 * @brief Encrypts a univariate integer polynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 *  TODO: implement
 *
 * @param module The underlying compute module
 * @param result The resulting GLWEGadget ciphertext
 * @param pk     A public key
 * @param m_univ The univariate plaintext
 *
 */
void glwegadget_public_encrypt(const MODULE* module, GLWEGadgetCiphertext* result, const GLWEPublicKey* pk,
                               const PolyUniv* m_univ);

/**
 * @brief Retreive a bivglwe from a GLWEGadget (halfGGSW for k = 1)
 *
 * @param glwegadget_ct     The GLWEGadget ciphertext
 * @param i                 The index of the GLWE to retrievie, from 1 to l_tilde included
 *
 */
VecBiv* glwegadget_extract_bivglwe(GLWEGadgetCiphertext* glwegadget_ct, uint64_t i);

/**
 * @brief Prepares a GLWEGadget for use in a half-external product
 *
 * @param module The underlying backend module
 * @param glwegadget_prep_ct The resulting prepared GLWEGadget
 * @param glwegad_ct The input GLWEGadget
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 *
 */
int glwegadget_prepare(const MODULE* module, GLWEGadgetCiphertextPrep* glwegadget_prep_ct,
                       const GLWEGadgetCiphertext* glwegad_ct);
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
