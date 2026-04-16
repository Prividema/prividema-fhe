#ifndef bivGGSW_CIPHERTEXT_H
#define bivGGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 *
 * In this header file, we define the structure representing bivariate GGSW ciphertext in both DFT and iDFT forms.
 */

#include "ggsw_params.h"
#include "spqlios_alias.h"

/**
 * @struct GGSWParams
 *
 * @brief GGSW Ciphertext in the coefficient space.
 */
typedef struct ggsw_ciphertext
{
	/// GGSW parameters.
	const GGSWParams* params;

	/// A matrix of size \f$n\_limbs\_tilde \times n\_limbs\f$ with coefficients that are in \f$\mathbb{Z}_n[X]\f$.
	MatBiv* mat;

} GGSWCiphertext;

/**
 * @brief Creates a Bivariate GGSW ciphertext filled with 0.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @retval • `NULL` if malloc failed inside the function.
 * @retval • The Allocated GGSW ciphertext filled with 0 otherwise.
 */
GGSWCiphertext* new_ggsw(const GGSWParams* params_ggsw);

/**
 * @brief Deletes a GGSW ciphertext.
 *
 * @param ggsw A Pointer to the GGSW parameters.
 */
void delete_ggsw(GGSWCiphertext* ggsw);

/**
 * @brief Get the Bivariate GLWE Ciphertext inside the GGSW with the given values.
 *
 * As a Bivariate GGSW's matrix is :
 *
 * - bivGLWE(DFT(-m * sk_j / 2^{kappa_tilde * i})) for j < k.
 * - bivGLWE(DFT(m / 2^{kappa_tilde * i}))         for j = k.
 *
 * The function takes j and i and get the associated bivGLWE.
 *
 * @param ggsw_ct    A Pointer to The GGSW ciphertext
 * @param j 		  The j-th component of the secret key.
 * @param i 		  The i in -m * sk_j / 2^{kappa_tilde * i}.
 *
 * @return A Pointer to the associated Bivariate GLWE.
 */
VecBiv* ggsw_retrieve_bivglwe(GGSWCiphertext* ggsw_ct, int64_t j, int64_t i);

/**
 * @brief Normalizes a GGSW ciphertext.
 *
 * @param module Additionnal information for backend.
 * @param result A Pointer to store the normalized GGSW ciphertext.
 * @param ggsw   A Pointer to the GGSW ciphertext to normalize.
 *
 * @retval -1 if an error occurs. In this case the error is from a malloc and perror is called.
 * @retval 0 otherwise.
 */
int normalize_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw);

/**
 * @brief Adds two GGSW ciphertexts with same parameters.
 *
 * @param result   A Pointer to store the result.
 * @param ggsw_lhs A Pointer to the left-hand side GGSW ciphertext.
 * @param ggsw_rhs A Pointer to the right-hand side GGSW ciphertext.
 */
void add_ggsw(GGSWCiphertext* result, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs);

/**
 * @brief Multiplies a GGSW ciphertext by a constant in \f$\mathbb{Z}_n[X]\f$.
 *
 * @param module    Additionnal information for backend.
 * @param result    A Pointer to store the result. Must be a GGSWCiphertext.
 * @param ggsw      A Pointer to the GGSW ciphertext.
 * @param cte       The constant in \f$\mathbb{Z}_n[X]\f$.
 *
 * @retval -1 if an error occurs. In this case the error is from a malloc and perror is called.
 * @retval 0 otherwise.
 */
int const_mult_ggsw(const MODULE* module, GGSWCiphertext* result, const GGSWCiphertext* ggsw, const PolyUnivDFT* cte);

// =============================================
// |                                           |
// |      			 DFT Form       		   |
// |                                           |
// | Same functions as above but in DFT Space. |
// |                                           |
// |     The parameters are the same but 	   |
// |   the ciphertexts are in the DFT space.   |
// =============================================

/**
 * @struct GGSWParams
 *
 * @brief GGSW Ciphertext in the DFT space.
 */
typedef struct ggsw_ciphertext_dft
{
	/// GGSW parameters.
	const GGSWParams* params;

	/// A matrix of size \f$n\_limbs\_tilde \times n\_limbs\f$ with coefficients that are in \f$\mathbb{C}_n[X]\f$.
	MatBivDFT* mat;

} GGSWCiphertextDFT;

/**
 * @brief Creates a Bivariate GGSW ciphertext in the DFT space filled with 0.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @retval • `NULL` if malloc failed inside the function.
 * @retval • The Allocated GGSW ciphertext filled with 0 otherwise.
 */
GGSWCiphertextDFT* new_ggsw_dft(const GGSWParams* params_ggsw);

/**
 * @brief Deletes a GGSW ciphertext in the DFT space.
 *
 * @param ggsw A Pointer to the GGSW parameters.
 */
void delete_ggsw_dft(GGSWCiphertextDFT* ggsw_dft);

/**
 * @brief Get the Bivariate GLWE Ciphertext inside the GGSW with the given values.
 *
 * As a Bivariate GGSW's matrix is :
 *
 * - bivGLWE(DFT(-m * sk_j / 2^{kappa_tilde * i})) if j < k.
 * - bivGLWE(DFT(m / 2^{kappa_tilde * i}))         if j = k.
 *
 * The function takes j and i and get the associated bivGLWE.
 *
 * @param ggsw_dft_ct    A Pointer to The GGSW ciphertext
 * @param j 		  The j-th component of the secret key.
 * @param i 		  The i in -m * sk_j / 2^{kappa_tilde * i}.
 *
 * @return A Pointer to the associated Bivariate GLWE in the DFT space.
 */
VecBivDFT* ggsw_retrieve_bivglwe_dft(GGSWCiphertextDFT* ggsw_dft_ct, int64_t j, int64_t i);

/**
 * @brief Adds two GGSW ciphertexts in the DFT space with same parameters.
 *
 * @param result   A Pointer to store the result.
 * @param ggsw_lhs A Pointer to the left-hand side GGSW ciphertext.
 * @param ggsw_rhs A Pointer to the right-hand side GGSW ciphertext.
 */
void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft,
                  const GGSWCiphertextDFT* ggsw_rhs_dft);

/**
 * @brief Multiplies a GGSW ciphertext in the DFT space by a constant in \f$\mathbb{Z}_n[X]\f$.
 *
 * @param module    Additionnal information for backend.
 * @param result    A Pointer to store the result. Must be a GGSWCiphertext.
 * @param ggsw      A Pointer to the GGSW ciphertext.
 * @param cte       The constant in \f$\mathbb{Z}_n[X]\f$. Must be in the DFT space.
 *
 * @retval -1 if an error occurs. In this case the error is from a malloc and perror is called.
 * @retval 0 otherwise.
 */
int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft);

#endif  // bivGGSW_CIPHERTEXT_H
