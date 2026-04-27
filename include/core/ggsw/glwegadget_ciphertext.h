#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

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
 * @brief Returns the number of coefficients in a GLWEGadgetCiphertext
 *
 * @param params_glwegadget The parameters
 *
 */
uint64_t glwegadget_coef_number(const GLWEGadgetParams* params_glwegadget);

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
#endif  // PARTIALGGSW_CIPHERTEXT_H
