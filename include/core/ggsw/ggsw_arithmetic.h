#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_ciphertext.h"
#include "glwe_ciphertext.h"
#include "glwe_key.h"
#include "glwegadget_key.h"
#include "spqlios_alias.h"

/**
 *
 * @file ggsw_arithmetic.h
 *
 * @brief This file contains arithmetic operations on GGSW ciphertexts
 *
 *
 */

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
 * @param module Additionnal information for backend.
 * @param result The GLWE result ciphertext.
 * @param glwe   The GLWE input ciphertext.
 * @param ggsw_prepared   The GGSW input ciphertext.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_external_product_to_dft(const MODULE* module, GLWECiphertextDFT* result, const GLWECiphertext* glwe,
                                 const GGSWCiphertextPrep* ggsw_prepared);

/**
 * @brief Computes the external product between a GLWE and a GGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The GLWE result ciphertext. Can be the same as the input
 * @param glwe   The GLWE input ciphertext.
 * @param ggsw_prepared   The GGSW input ciphertext, in prepared format.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_external_product(const MODULE* module, GLWECiphertext* result, const GLWECiphertext* glwe,
                          const GGSWCiphertextPrep* ggsw_prepared);

/**
 * @brief Expands a packed GLWEGadget into a set of corresponding GGSWs
 *
 *
 * @param module              The backend module
 * @param results             An array with pointers to the res_size output GGSWs
 * @param res_size            The number of (non-zero) coefficients in the packed GLWEGadget. Equivalently, the number of output GGSWs.
 * @param l_tilde             The l_tilde with which the packed GLWEGadget was packed
 * @param packed_glwegadget   The packed GLWEGadget
 * @param auto_key_collection The key collection for trace expansion automorphisms
 * @param sk_encryptions      An array of pointers to encryptions of GGSW(-sk_i) for i=1..k, used to convert GLWEGadgets into GGSWs
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int packed_glwegadget_trace_expand_ggsw(const MODULE* module, GGSWCiphertext** results, int res_size, int l_tilde,
                                        const GLWECiphertext* packed_glwegadget,
                                        const GLWEAutomorphismKeyCollection* auto_ksks,
                                        const GGSWCiphertextPrep** sk_encryptions);

/**
 * @brief Expands a packed GLWEGadget into a set of corresponding GGSWs
 *
 * See https://github.com/Prividema/prividema-fhe/pull/64 for an explanation on what a packed GLWEGadget is.
 * This variant expands the GLWEGadget and uses the GGSWCiphertextPrep encryptions of -sk_i to convert it
 * into a GGSW.
 *
 * @param module              The backend module
 * @param results              An array with pointers to the res_size output prepared GGSWs.
 *                             If NULL pointers are provided, they will be created and allocated by this function, and responsibility
 *                             tranferred to the caller. This is the recommended mode of operation
 * @param res_size             The number of (non-zero) coefficients in the packed GLWEGadget. Equivalently, the number of output GGSWs.
 * @param l_tilde              The l_tilde with which the packed GLWEGadget was packed
 * @param packed_glwegadget    The packed GLWEGadget
 * @param auto_key_collection  The KSK collection for trace expansion automorphisms
 * @param sk_encryptions       An array of pointers to encryptions to GGSW(-sk_i) for i=1..k, used to convert GLWEGadgets into GGSWs
 *
 *
 * @note It is recommended to call this function without having pre-allocated the results and therefore letting this function perform
 *       their allocation. Due to the algorithm used for trace expansion, it is required to expand an entire trace at once.
 *       However, that is done in the unprepared space.
 *       If this function is called with preallocated results, at some point, memory will be allocated for the entire unprepared
 *       trace as well as the entire prepared trace. This results in almost double the peak memory consumption that is really needed.
 *       If the function is allowed to allocate the output, it will do so one element of the trace at a time, while at the same time
 *       it will deallocate the elements of the unprepared trace that it no longer requires. Therefore, the peak memory
 *       consumption will be the size of a trace plus one element (the one being prepared) instead of two traces (one prepared, one unprepared).
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int packed_glwegadget_trace_expand_ggsw_prepared(const MODULE* module, GGSWCiphertextPrep** results,
                                                 const GGSWParams* params_ggsw, int res_size, int l_tilde,
                                                 const GLWECiphertext* packed_glwegadget,
                                                 const GLWEAutomorphismKeyCollection* auto_key_collection,
                                                 const GGSWCiphertextPrep** sk_encryptions);
#endif  // bivGGSW_H
