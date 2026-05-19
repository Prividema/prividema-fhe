#ifndef bivGGSW_CIPHERTEXT_H
#define bivGGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 *
 * In this header file, we define the structure representing bivariate GGSW ciphertext in both DFT and coefficient forms.
 */

#include "ggsw_params.h"
#include "glwe_key.h"
#include "maths_structures.h"
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
 * @param ggsw_ct     A Pointer to The GGSW ciphertext
 * @param sk_idx      The j-th component of the secret key.
 * @param prec_lvl    The i in \f$ -m · sk_j · 2^{-i\kappatilde} \f$
 *
 * @return A Pointer to the associated Bivariate GLWE.
 */
VecBiv* ggsw_retrieve_bivglwe(GGSWCiphertext* ggsw_ct, int64_t sk_idx, int64_t prec_lvl);

/**
 * @brief Encrypts the \ZnX message m into the bivGGSW ciphertext res according to the
 * parameters in the result object.
 *
 * @param module      Additionnal information for backend.
 * @param result      The encrypted message.
 * @param sk_prep      The prepared secret key.
 * @param m_univ      The univariate message (in coefficient space).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int ggsw_secret_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWESecretKeyPrepared* sk_prep,
                        const PolyUniv* m_univ);

/**
 * @brief Encrypts univariate \ZnX message m_univ into GGSW ciphertext res with parameters params
 *
 * TODO: implement in a future release
 *
 * @param module The backend module
 * @param result The resulting GGSWCiphertext
 * @param pk     The public key used for encryption
 * @param m_univ The univariate message to encrypt (in coefficient space)
 *
 * */
int ggsw_public_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWEPublicKey* pk, const PolyUniv* m_univ);

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

#endif  // bivGGSW_CIPHERTEXT_H
