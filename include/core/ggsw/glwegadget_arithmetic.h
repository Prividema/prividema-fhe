#ifndef PARTIAL_GGSW_H
#define PARTIAL_GGSW_H

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwegadget_ciphertext.h"
#include "glwegadget_key.h"

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

/**
 *
 * @brief Creates a key-switching-key (KSK) for an automorphism of degree automorphism_p from
 * the provided secret key
 *
 * @param module The backend module
 * @param automorphism_ksk The output generated key-switching-key
 * @param glwe_key The input prepared secret key
 * @param automorphism_p The degree of the automorphism. Can be positive or negative,
 * the automorphism is only well-defined if p is odd
 *
 *
 */
int prepare_automorphism_key(const MODULE* module, GLWEAutomorphismKSK* automorphism_ksk,
                             const GLWESecretKeyPrepared* glwe_key, int automorphism_p);

/**
 * @brief Performs the automorphism of a GLWE ciphertext , ie, sets its polynomials a(X, Y) to a(X^p, Y)
 *
 * @param module             The backend module
 * @param result             The resulting GLWE ciphertext
 * @param automorphism_ksk An encryption of the secret key after having applied the automorphism to it
 * @param glwe The input ciphertext
 * @param automorphism_p The p value for the automorphism
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_automorphism(const MODULE* module, GLWECiphertext* result, const GLWEAutomorphismKSK* automorphism_ksk,
                            const GLWECiphertext* glwe, int automorphism_p);

int glwe_trace_expand(const MODULE* module, GLWECiphertext** results, int res_size, const GLWECiphertext* glwe_ct,
                      const GLWEAutomorphismKSKCollection* ksks);

int packed_glwegadget_trace_expand(const MODULE* module, GLWEGadgetCiphertext** results, int res_size, int l_tilde,
                                   const GLWECiphertext* packed_glwegadget,
                                   const GLWEAutomorphismKSKCollection* auto_ksks);
/**
 *
 * The function can allocate results by itself, to allow for reduced peak memory consumption (since trace expansion can be demanding)
 * In order to do that, simply have results be an array of null pointers
 *
 *
 */
int packed_glwegadget_trace_expand_prepared(const MODULE* module, GLWEGadgetCiphertextPrep** results, int res_size,
                                            const GLWECiphertext* packed_glwegadget,
                                            const GLWEAutomorphismKSKCollection* auto_ksks);
#endif  // !DEBUG
