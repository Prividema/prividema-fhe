#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include "bivariate_polynomial.h"
#include "ggsw_params.h"

typedef struct partialggsw_ciphertext
{
	const GLWEGadParams* params;
	MatBiv* mat;  // matrix of size n_limbs x l_tilde
} GLWEGadCiphertext;

/**
 * @brief creates a new GLWEGadget ciphertext
 *
 * @params params The parameters
 *
 */
GLWEGadCiphertext* new_glwegad(const GLWEGadParams* params);

/**
 * @brief Returns the number of coefficients in a GLWEGadgetCiphertext
 *
 * @params params_glwegad The parameters
 *
 */
size_t glwegad_coef_number(const GLWEGadParams* params_glwegad);

/**
 * @brief Deletes a GLWEGadget ciphertext
 *
 * Remember that it does not own the underlying parameters
 *
 * @params glwegad_ct The ciphertext to delete
 *
 */
void delete_glwegad(GLWEGadCiphertext* glwegad_ct);

typedef struct partialggsw_ciphertext_prepared
{
	const GLWEGadParams* params;
	MatBivDFT* mat;  // matrix of size n_limbs x l_tilde
} GLWEGadCiphertextPrep;

/**
 * @brief creates a new partialggsw ciphertext
 *
 * @params params The parameters
 *
 */
GLWEGadCiphertextPrep* new_glwegad_prep(const GLWEGadParams* params);

/**
 * @brief Deletes a PartialGGSWCiphertext
 *
 * Remember that it does not own the underlying parameters
 *
 * @params glwegad_ct The ciphertext to delete
 *
 */
void delete_glwegad_prep(GLWEGadCiphertext* glwegad_prep_ct);

#endif  // PARTIALGGSW_CIPHERTEXT_H
