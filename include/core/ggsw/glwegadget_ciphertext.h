#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include "bivariate_polynomial.h"
#include "ggsw_params.h"

typedef struct partialggsw_ciphertext
{
	const GLWEGadgetParams* params;
	MatBiv* mat;  // matrix of size n_limbs x l_tilde
} GLWEGadgetCiphertext;

/**
 * @brief creates a new GLWEGadget ciphertext
 *
 * @params params The parameters
 *
 */
GLWEGadgetCiphertext* new_glwegadget(const GLWEGadgetParams* params);

/**
 * @brief Returns the number of coefficients in a GLWEGadgetCiphertext
 *
 * @params params_glwegadget The parameters
 *
 */
size_t glwegadget_coef_number(const GLWEGadgetParams* params_glwegadget);

/**
 * @brief Deletes a GLWEGadget ciphertext
 *
 * Remember that it does not own the underlying parameters
 *
 * @params glwegadget_ct The ciphertext to delete
 *
 */
void delete_glwegadget(GLWEGadgetCiphertext* glwegadget_ct);

typedef struct partialggsw_ciphertext_prepared
{
	const GLWEGadgetParams* params;
	MatBivDFT* mat;  // matrix of size n_limbs x l_tilde
} GLWEGadgetCiphertextPrep;

/**
 * @brief creates a new partialggsw ciphertext
 *
 * @params params The parameters
 *
 */
GLWEGadgetCiphertextPrep* new_glwegadget_prep(const GLWEGadgetParams* params);

/**
 * @brief Deletes a PartialGGSWCiphertext
 *
 * Remember that it does not own the underlying parameters
 *
 * @params glwegadget_prep_ct The ciphertext to delete
 *
 */
void delete_glwegadget_prep(GLWEGadgetCiphertext* glwegadget_prep_ct);

#endif  // PARTIALGGSW_CIPHERTEXT_H
