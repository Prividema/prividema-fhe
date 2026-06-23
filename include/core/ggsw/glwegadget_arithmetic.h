#ifndef PARTIAL_GGSW_H
#define PARTIAL_GGSW_H

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"

/**
 * @file glwegadget_arithmetic.h
 *
 * @brief Arithmetic operations for GLWEGadgets
 *
 */

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

int glwegadget_half_prod_prepared_to_dft(const MODULE* module, GLWECiphertextDFT* result_dft,
                                         const GLWEGadgetCiphertextPrep* glwegadget_prep_ct, const PolyBivPrep* a_prep);
/**
 * @brief Creates a key-switching-key (KSK) for an automorphism of degree automorphism_p from
 * the provided secret key
 *
 * @param module              The backend module
 * @param automorphism_ksk    The output generated key-switching-key
 * @param glwe_key            The input prepared secret key
 * @param automorphism_p      The degree of the automorphism. Can be positive or negative,
 *                            the automorphism is only well-defined if p is odd
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int prepare_automorphism_key(const MODULE* module, GLWEAutomorphismKSK* automorphism_ksk,
                             const GLWESecretKeyPrepared* glwe_key, int automorphism_p);
/**
 * @brief Prepares a KSK for GLWE keyswitching from old_key to new_key
 *
 * @param module        The backend module
 * @param ksk           The output KSK
 * @param old_key       The old key being switched from
 * @param new_key       The new key being switched to
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int prepare_ksk(const MODULE* module, GLWEAutomorphismKSK* ksk, const GLWESecretKeyPrepared* new_key,
                const GLWESecretKeyPrepared* old_key);
/**
 * @brief Performs the automorphism of a GLWE ciphertext, i.e. sets its polynomials a(X, Y) to a(X^p, Y)
 *
 * @param module             The backend module
 * @param result             The resulting GLWE ciphertext
 * @param automorphism_ksk   An encryption of the secret key after having applied the automorphism to it.
 *                           The function retrieves the p value from it.
 * @param glwe               The input ciphertext
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_automorphism(const MODULE* module, GLWECiphertext* result, const GLWEAutomorphismKSK* automorphism_ksk,
                            const GLWECiphertext* glwe);
/**
 * @brief Performs a key-switch on a GLWE using an already prepared KSK
 *
 * @param module             The backend module
 * @param result             The resulting GLWE keyed with the new key
 * @param automorphism_ksk   A prepared KSK from the old key to the new key
 * @param glwe_ct            The input GLWE keyed with the old key
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwe_to_glwe_keyswitch(const MODULE* module, GLWECiphertext* result, const GLWEAutomorphismKSK* automorphism_ksk,
                           const GLWECiphertext* glwe_ct);

/**
 * @brief Expands a GLWE with at most res_size non-zero coefficents into a res_size GLWEs encrypting each coefficient
 *
 * More specifically, expands GLWE(a_0 + a_1 * x + ... + a_d * x^d) into GLWE(m * a_0) ... GLWE(m * a_d) with d = res_size - 1
 * and m = ceil(log2(res_size))
 *
 * The algorithm used as well as a better explanation can be found in
 * https://github.com/Prividema/prividema-fhe/pull/64#issuecomment-4648948222
 *
 *
 * @param module The backend module
 * @param results An array with pointers to the res_size output GLWEs.
 * @param res_size The number of (non-zero) coefficients in the GLWE.
 * @param glwe_ct The packed GLWEGadget
 * @param ksks The KSK collection for trace expansion automorphisms
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwe_trace_expand(const MODULE* module, GLWECiphertext** results, int res_size, const GLWECiphertext* glwe_ct,
                      const GLWEAutomorphismKSKCollection* ksks);

/**
 * @brief Expands a packed GLWEGadget into a set of corresponding GLWEGadgets
 *
 * See https://github.com/Prividema/prividema-fhe/pull/64 for an explanation on "packed" GLWEGadgets
 *
 * @param module The backend module
 * @param results An array with pointers to the res_size output GLWEGadgets.
 * @param res_size The number of (non-zero) coefficients in the packed GLWEGadget. Equivalently, the number of output GLWEGadgets
 * @param l_tilde The l_tilde value used when packing and encrypting the GLWEGadget
 * @param packed_glwegadget The packed GLWEGadget
 * @param auto_ksks The KSK collection for trace expansion automorphisms
 *
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int packed_glwegadget_trace_expand(const MODULE* module, GLWEGadgetCiphertext** results, int res_size, int l_tilde,
                                   const GLWECiphertext* packed_glwegadget,
                                   const GLWEAutomorphismKSKCollection* auto_ksks);
/**
 * @brief Expands a packed GLWEGadget into a set of corresponding GLWEGadgets
 *
 * See https://github.com/Prividema/prividema-fhe/pull/64 for an explanation on "packed" GLWEGadgets
 *
 * @param module The backend module
 * @param results An array with pointers to the res_size output prepared GLWEGadgets
 *                If NULL pointers are provided, they will be created and allocated by this function, and responsibility
 *                tranferred to the caller. This is the recommended mode of operation
 * @param res_size The number of (non-zero) coefficients in the packed GLWEGadget. Equivalently, the number of output GLWEGadgets
 * @param packed_glwegadget The packed GLWEGadget
 * @param auto_ksks The KSK collection for trace expansion automorphisms
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
int packed_glwegadget_trace_expand_prepared(const MODULE* module, GLWEGadgetCiphertextPrep** results, int res_size,
                                            int l_tilde, const GLWECiphertext* packed_glwegadget,
                                            const GLWEAutomorphismKSKCollection* auto_ksks);

int packed_glwegadget_trace_expand_prepared_single(const MODULE* module, GLWEGadgetCiphertextPrep* results,
                                                   const GLWEGadgetParams* params_glwegad, int res_size, int l_tilde,
                                                   const GLWECiphertext* packed_glwegadget,
                                                   const GLWEAutomorphismKSKCollection* auto_ksks);
#endif  // !DEBUG
