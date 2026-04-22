#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include "bivariate_polynomial.h"
#include "ggsw_params.h"

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

#endif  // PARTIALGGSW_CIPHERTEXT_H
