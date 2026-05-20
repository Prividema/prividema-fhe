#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwegadget_key.h"
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
 * @param module    Additional information for backend.
 * @param result    The resulting GGSWCiphertext
 * @param ggsw_lhs  One of the input GGSWCiphertext
 * @param ggsw_rhs  The other input GGSWCiphertext
 */
void add_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw_lhs,
              const GGSWCiphertext* ggsw_rhs);

/**
 * @brief Multiplies a GGSW ciphertext by a constant in \ZnX
 *
 * @param module    Additional information for backend.
 * @param result    The resulting GGSWCiphertext
 * @param ggsw      The input GGSWCiphertext.
 * @param cte       The constant in \ZnX
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* cte);

/**
 * @brief Computes the external product between a GLWE and a GGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The GLWE result ciphertext. Can be the same as the input
 * @param glwe   The GLWE input ciphertext.
 * @param ggsw   The GGSW input ciphertext.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_unprepared_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                                     const GGSWCiphertext* ggsw);

/**
 * @brief Computes the external product between a GLWE and a GGSW.
 *
 * TODO: change input to DFT form one vmp dft_to_dft is implemented
 *
 * @param module Additionnal information for backend.
 * @param result The GLWE result ciphertext.
 * @param glwe   The GLWE input ciphertext.
 * @param ggsw   The GGSW input ciphertext.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_external_product_to_dft(const MODULE* module, GLWECiphertextDFT* result, const GLWECiphertext* glwe,
                                 const GGSWCiphertextPrep* ggsw_prepared);

int ggsw_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                          const GGSWCiphertextPrep* ggsw_prepared);

/**
 *
 * Expands a packed GLWEGadget query into a collection of GGSWs
 *
 *
 *
 */
int packed_glwegadget_trace_expand_ggsw(const MODULE* module, GGSWCiphertext** results, int res_size, int l_tilde,
                                        const GLWECiphertext* packed_glwegadget,
                                        const GLWEAutomorphismKSKCollection* auto_ksks,
                                        const GGSWCiphertextPrep** sk_encryptions);

#endif  // bivGGSW_H
