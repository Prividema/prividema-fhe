#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include <stdint.h>

#include "bivariate_polynomial.h"
#include "ggsw_params.h"
#include "glwe_key.h"

/**
 * @brief A GLWEGadget ciphertext
 */
typedef struct glwegadget_ciphertext
{
	const GLWEGadgetParams* params;
	MatBiv* mat;  ///< matrix of size l_tilde * n_limbs, each row the contents of a GLWECiphertext
} GLWEGadgetCiphertext;

/**
 * @brief creates a new GLWEGadgetCiphertext
 *
 * @param params The parameters
 *
 */
GLWEGadgetCiphertext* new_glwegadget(const GLWEGadgetParams* params);

/**
 * @brief Deletes a GLWEGadget ciphertext
 *
 * Remember that it does not own the underlying parameters
 *
 * @param glwegadget_ct The ciphertext to delete
 *
 */
void delete_glwegadget(GLWEGadgetCiphertext* glwegadget_ct);

/**
 * @brief A GLWEGadget ciphertext that has been preprocessed for use in a half-external product
 */
typedef struct glwegadget_ciphertext_prepared
{
	const GLWEGadgetParams* params;
	MatBivDFT*
	    mat;  ///< Prepared GLWEGadgetCiphertext matrix. Its layout is an internal implementation detail of the underlying compute engine
} GLWEGadgetCiphertextPrep;

/**
 * @brief creates a new prepared GLWEGadget
 *
 * @param params The parameters
 *
 */
GLWEGadgetCiphertextPrep* new_glwegadget_prep(const GLWEGadgetParams* params);

/**
 * @brief Deletes a prepared GLWEGadget
 *
 * Remember that it does not own the underlying parameters
 *
 * @param glwegadget_prep_ct The ciphertext to delete
 *
 */
void delete_glwegadget_prep(GLWEGadgetCiphertextPrep* glwegadget_prep_ct);

/**
 * @brief Encrypts a univariate integer polynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 * @param module The underlying compute module
 * @param result The resulting GLWEGadget ciphertext
 * @param sk_prep A prepared secret key
 * @param m_univ The univariate plaintext
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwegadget_secret_encrypt(const MODULE* module, GLWEGadgetCiphertext* result, const GLWESecretKeyPrepared* sk_prep,
                              const PolyUniv* m_univ);

/**
 * @brief Encrypts a univariate integer polynomial into a GLWEGadget
 * (also known as halfGGSW in the k=1 case)
 *
 *  TODO: not yet implemented
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
 *
 * Encrypts a univariate polynomial with at most d non-zero coefficients
 * into a "packed GLWEGadget" GLWE ciphertext (see @ref )
 *
 *
 *
 */
int glwegadget_packed_secret_encrypt(const MODULE* module, GLWECiphertext* result,
                                     const GLWEGadgetParams* params_glwegad, const GLWESecretKeyPrepared* sk_prep,
                                     const PolyUniv* m_univ, uint64_t d);

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

#endif  // PARTIALGGSW_CIPHERTEXT_H
