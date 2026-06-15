#ifndef bivGGSW_CIPHERTEXT_H
#define bivGGSW_CIPHERTEXT_H

/**
 * @file ggsw_ciphertext.h
 *
 * @brief Header containing GGSW ciphertext definition, allocation and deallocation as well as encryption
 *
 * Decryption is not implemented nor is is planned.
 * See the comments in glwegadget_ciphertext.h for the reasoning, which is the same for GGSWs.
 *
 *
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
	/// The GGSW parameters
	const GGSWParams* params;

	/// A flattened (in row-major fashion) matrix of size
	/// n_limbs_tilde \f$ \cdot \f$ n_limbs with each element a polynomial \ZnX,
	/// representing a GGSW ciphertext (as a series of GLWE ciphertexts).
	/// It is NOT a concatenation of GLWEGadgets, as the strided
	/// representation of them that this structure uses makes it have a prefix
	/// property on the rows as well as the columns.
	///
	/// See \ref ggsw_encoding "the section on GGSW ciphertext encoding" for a
	/// more in-depth explanation.
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
 * @brief Get the bivariate GLWE Ciphertext inside the GGSW with the given values.
 *
 * A bivariate GGSW's matrix is :
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

PolyBiv ggsw_flattened_biv(const GGSWCiphertext* ggsw_ct);

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
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 * */
int ggsw_public_encrypt(const MODULE* module, GGSWCiphertext* result, const GLWEPublicKey* pk, const PolyUniv* m_univ);

// =============================================
// |                                           |
// |      			 Prepared Form       		   |
// |                                           |
// =============================================

/**
 * @brief GGSW Ciphertext in the DFT space.
 */
typedef struct ggsw_ciphertext_prepared
{
	/// GGSW parameters.
	const GGSWParams* params;

	/// Prepared matrix (opaque representation)
	MatBivDFT* mat;

} GGSWCiphertextPrep;

/**
 * @brief Creates a zero-initialized GGSWCiphertextPrep
 *
 * @param params_ggsw The GGSW parameters.
 *
 * @retval NULL if malloc failed inside the function.
 * @retval non-NULL The Allocated GGSW ciphertext filled with 0 otherwise.
 */
GGSWCiphertextPrep* new_ggsw_prep(const GGSWParams* params_ggsw);

/**
 * @brief Deletes a prepared GGSW ciphertext
 *
 * @param ggsw_dft A Pointer to the prepared GGSW
 */
void delete_ggsw_prep(GGSWCiphertextPrep* ggsw_dft);

/**
 * @brief Prepares a GGSW ciphertext for an efficient external product
 *
 * @param module The backend module
 * @param ggsw_prepared The output prepared GGSW ciphertext
 * @param ggsw_ct The input unprepared GGSW ciphertext
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int ggsw_prepare(const MODULE* module, GGSWCiphertextPrep* ggsw_prepared, const GGSWCiphertext* ggsw_ct);

#endif  // bivGGSW_CIPHERTEXT_H
