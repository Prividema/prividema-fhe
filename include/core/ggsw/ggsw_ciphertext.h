#ifndef bivGGSW_CIPHERTEXT_H
#define bivGGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 *
 * In this header file, we define the structure representing bivariate GGSW ciphertext in both DFT and coefficient forms.
 */

#include "ggsw_params.h"
#include "spqlios_alias.h"

/**
 * @brief GGSW Ciphertext in the coefficient space.
 */
typedef struct ggsw_ciphertext
{
	/// GGSW parameters.
	const GGSWParams* params;

	/// A matrix of size n_limbs_tilde \f$·\f$ n_limbs with coefficients that are in \ZnX
	MatBiv* mat;

} GGSWCiphertext;

/**
 * @brief Creates a Bivariate GGSW ciphertext filled with 0.
 *
 * @param params_ggsw A Pointer to the GGSW parameters.
 *
 * @retval NULL if malloc failed inside the function.
 * @retval non-NULL The Allocated GGSW ciphertext filled with 0 otherwise.
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
 * - bivGLWE(DFT(-m * sk_j / 2^{\kappa_tilde * i})) for j < k.
 * - bivGLWE(DFT(m / 2^{\kappa_tilde * i}))         for j = k.
 *
 * The function takes j and i and get the associated bivGLWE.
 *
 * @param ggsw_ct    A Pointer to The GGSW ciphertext
 * @param j 		  The j-th component of the secret key.
 * @param i 		  The i in \f$ -m · sk_j · 2^{-i\kappatilde} \f$
 *
 * @return A Pointer to the associated Bivariate GLWE.
 */
VecBiv* ggsw_retrieve_bivglwe(GGSWCiphertext* ggsw_ct, int64_t j, int64_t i);

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
 * @param result    The resulting GGSWCiphertext
 * @param ggsw_lhs  One of the input GGSWCiphertext
 * @param ggsw_rhs  The other input GGSWCiphertext
 */
void add_ggsw(GGSWCiphertext* result, const GGSWCiphertext* ggsw_lhs, const GGSWCiphertext* ggsw_rhs);

/**
 * @brief Multiplies a GGSW ciphertext by a constant in \ZnX
 *
 * @param module    Additionnal information for backend.
 * @param result    The resulting GGSWCiphertext
 * @param ggsw      The input GGSWCiphertext.
 * @param cte       The constant in \ZnX
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
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
 * @brief GGSW Ciphertext in the DFT space.
 */
typedef struct ggsw_ciphertext_dft
{
	/// GGSW parameters.
	const GGSWParams* params;

	/// A matrix of polynomials in the DFT domain (opaque representation)
	MatBivDFT* mat;

} GGSWCiphertextDFT;

/**
 * @brief Creates a zero-initialized GGSWCiphertextDFT.
 *
 * @param params_ggsw The GGSW parameters.
 *
 * @retval NULL if malloc failed inside the function.
 * @retval non-NULL The Allocated GGSW ciphertext filled with 0 otherwise.
 */
GGSWCiphertextDFT* new_ggsw_dft(const GGSWParams* params_ggsw);

/**
 * @brief Deletes a GGSW ciphertext in the DFT space.
 *
 * @param ggsw_dft A Pointer to the DFT GGSW
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
 * @param ggsw_dft_ct   The GGSW ciphertext
 * @param j 		  The j-th component of the secret key.
 * @param i 		  The i in \f$ -m · sk_j · 2^{-i\kappatilde} \f$
 *
 * @return A pointer to the associated bivGLWEDFT
 */
VecBivDFT* ggsw_retrieve_bivglwe_dft(GGSWCiphertextDFT* ggsw_dft_ct, int64_t j, int64_t i);

/**
 * @brief Adds two GGSW ciphertexts in the DFT space  that have the same parameters
 *
 * @param result_dft   The resulting GGSWCiphertextDFT
 * @param ggsw_lhs_dft One of the input GGSWCiphertextDFT.
 * @param ggsw_rhs_dft The other GGSWCiphertextDFT.
 */
void add_ggsw_dft(GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_lhs_dft,
                  const GGSWCiphertextDFT* ggsw_rhs_dft);

/**
 * @brief Multiplies a GGSW ciphertext in the DFT space by a polynomial from \ZnX also in DFT space
 *
 * @param module    Additionnal information for backend.
 * @param result_dft The output GGSWCiphertextDFT
 * @param ggsw_dft      The input GGSWCiphertextDFT.
 * @param u_dft       The \ZnX polynomial in the DFT space.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int const_mult_ggsw_dft(const MODULE* module, GGSWCiphertextDFT* result_dft, const GGSWCiphertextDFT* ggsw_dft,
                        const PolyUnivDFT* u_dft);

#endif  // bivGGSW_CIPHERTEXT_H
