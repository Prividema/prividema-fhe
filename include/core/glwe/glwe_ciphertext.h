#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>

#include "backend_types.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "maths_structures.h"
#include "bivariate_polynomial.h"

/**
 * @file glwe_ciphertext.h
 *
 * In this header file, we define the structures representing bivariate GLWE ciphertext in both DFT and coefficient forms,
 * encryption operations to them, decryption, allocation and deallocation and accessors to their components.
 */

// =============================================
// |                                           |
// |             GLWE Ciphertext               |
// |                                           |
// =============================================

/**
 * @brief A GLWE ciphertext, encoded in bivariate base-2K form
 */
typedef struct glwe_ciphertext
{
	const GLWEParams* params;  ///< GLWE parameters

	/**
	 * @brief Represents a vector of size nb_limbs with coefficients that are in \ZnX
	 *
	 * Data layout is limb-major, ie, all the coefficients for a power of \f$ 2^{iK} \f$
	 *  are stored contiguously after the ones from \f$ 2^{(i-1)K} \f$ and before \f$ 2^{(i+1)K} \f$
	 * (see dedicated document/page for details)
	 */
	int64_t* vec;
} GLWECiphertext;

/**
 * @brief Creates a GLWE ciphertext, filled with 0. This function allocates the memory.
 * 
 * @param params_glwe The GLWE parameters.
 * 
 * @return The newly allocated object
 */
GLWECiphertext* new_glwe(const GLWEParams* params_glwe);

/**
 * @brief Deletes a GLWE ciphertext. This function frees the memory.
 * 
 * @param glwe The GLWE ciphertext to delete. 
 */
void delete_glwe(GLWECiphertext* glwe);

/**
 * @brief Copies a GLWE ciphertext.
 * 
 * This copy function truncates if dst is smaller than src.
 * This behavior is completely intentional as it means that we drop the less significant digits.
 * If dst is larger than src, the extra (less significant) limbs of dst are set to 0.
 *
 * @param dst The destination GLWE ciphertext.
 * @param src The source GLWE ciphertext.
 */
void glwe_copy(GLWECiphertext* dst, const GLWECiphertext* src);

/**
 * @brief Returns a PolyBiv view of the pos'th element (in a_0, ..., a_k = b) in a GLWE
 *
 * @param glwe_ct Source GLWE
 * @param pos Position to create a view of
 *
 * @returns A PolyBiv view (does NOT own memory) of the pos'th element in the GLWE
 */
PolyBiv glwe_extract_poly_view(const GLWECiphertext* glwe_ct, uint64_t pos);

/**
 * @brief Flat Znx view of the whole ciphertext (all limbs of all elements, in memory order), for element-wise operations.
 *
 * @param glwe_ct Source GLWE
 *
 * @returns A ZnxVec view (does NOT own memory) of all the limbs of the GLWE, as a flat vector of polynomials
 */
ZnxVec glwe_as_znx_vec(const GLWECiphertext* glwe_ct);

/**
 * @brief Encrypts a phase (message + noise) and puts it in result.
 *
 * @param module Additional information for backend.
 * @param res The result bivariate ciphertext.
 * @param sk_prep The prepared secret key
 * @param phase message + noise.
 *
 * @retval -1 if an error occurs.
 * @retval  0 otherwise.
 */
int glwe_secret_encrypt_phase(const PvdaBackend* module, GLWECiphertext* res, const GLWESecretKeyPrepared* sk_prep,
                              const PolyBiv* phase);

/**
 * @brief Encrypts a univariate message in RnX form
 *
 * @param module Additional information for backend.
 * @param res The result bivariate ciphertext.
 * @param sk_prep The prepared secret key
 * @param m_univ_rnx  The univariate secret message to encrypt
 *
 * @retval -1 if an error occurs.
 * @retval  0 otherwise.
 */
int glwe_secret_encrypt_rnx(const PvdaBackend* module, GLWECiphertext* res, const GLWESecretKeyPrepared* sk_prep,
                            const PolyUnivRnX* m_univ_rnx);

/**
 * @brief Encrypts a univariate message in TnX form
 *
 * @param module Additional information for backend.
 * @param res The result bivariate ciphertext.
 * @param sk_prep The prepared secret key
 * @param m_univ_tnx  The univariate secret message to encrypt
 *
 * @retval -1 if an error occurs.
 * @retval  0 otherwise.
 */
int glwe_secret_encrypt_tnx(const PvdaBackend* module, GLWECiphertext* res, const GLWESecretKeyPrepared* sk_prep,
                            const PolyUnivTnX* m_univ_tnx);

/**
 * @brief Decrypts the ciphertext into the
 * bivariate phase (message + noise) and puts it in result.
 *
 * @param module Additional information for backend.
 * @param res The bivariate phase.
 * @param sk_prep The prepared secret key
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval -1 if an error occurs.
 * @retval  0 otherwise.
 */
int glwe_secret_decrypt(const PvdaBackend* module, PolyBiv* res, const GLWESecretKeyPrepared* sk_prep,
                        const GLWECiphertext* glwe);

// =============================================
// |                                           |
// |           GLWE Ciphertext DFT             |
// |                                           |
// =============================================

/**
 * @brief A GLWE ciphertext in the DFT domain, encoded in bivariate base-2K form.
 */
typedef struct glwe_ciphertext_dft
{
	const PvdaBackend* module;	///< Backend that owns vec
	const GLWEParams* params;	///< GLWE parameters
	BE_PolyVecDFT* vec;  		///< Opaque DFT-domain vector, one entry per limb of the GLWE.
} GLWECiphertextDFT;

/**
 * @brief Creates a new empty GLWE ciphertext in DFT domain. This function allocates the memory.
 * @param params_glwe The GLWE parameters.
 * @return A pointer to the new object or NULL in case of a failure.
 */
GLWECiphertextDFT* new_glwe_dft(const PvdaBackend* module, const GLWEParams* params_glwe);

/**
 * @brief Deletes a GLWE ciphertext, but not its parameters. This function frees the memory.
 * @param glwe The GLWE ciphertext to delete.
 */
void delete_glwe_dft(GLWECiphertextDFT* glwe);

// =============================================
// |                                           |
// |        GLWE Ciphertext Transforms         |
// |                                           |
// =============================================

/**
 * @brief Convert a GLWECiphertext into its DFT representation.
 *
 * @param module  The backend module.
 * @param res_dft The resulting GLWE ciphertext in the DFT domain.
 * @param glwe_ct The input GLWE ciphertext in coeff space.
 *
 * @retval -1 If an error occurs.
 * @retval 0 otherwise.
 */
int glwe_coef_to_dft(const PvdaBackend* module, GLWECiphertextDFT* res_dft, const GLWECiphertext* glwe_ct);

/**
 *
 * @brief Convert a GLWECiphertextDFT into its coefficient (non-DFT) representation.
 *
 * @param module   The backend module.
 * @param res_ct   The resulting GLWE ciphertext in coeff space.
 * @param glwe_dft The input GLWE ciphertext in the DFT domain.
 *
 * @retval -1 If an error occurs.
 * @retval 0 otherwise.
 */
int glwe_dft_to_coef(const PvdaBackend* module, GLWECiphertext* res_ct, const GLWECiphertextDFT* glwe_dft);

// =============================================
// |                                           |
// |        	  Debug Functions              |
// |                                           |
// =============================================

/**
 * @brief Print function for debugging.
 *
 * Intended to inspect GLWE content for debugging.
 * It takes as input a secret key to decrypt. For debugging, a global variable could
 * be used to store the reference.
 *
 * Prints the n first (lowest degree) coefficients of the decryption.
 * Additionally, if shft is non-zero, the resulting coefficients will be bit-shifted
 * by that amount, which is intended to be used with BFV-like encodings.
 *
 * Only able to output the 64 MSBs of the encrypted message/phase
 *
 * @param module The backend module
 * @param glwe The glwe to print
 * @param sk_prep The prepared secret key for the glwe
 * @param n How many coefficients to print
 * @param shft How much to right-shift the end result to truncate noise in BFV-like encodings
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise
 */
int glwe_print_coefs(const PvdaBackend* module, const GLWECiphertext* glwe, const GLWESecretKeyPrepared* sk_prep, int n, int shft);

#endif  // GLWE_CIPHERTEXT_H
