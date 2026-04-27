#ifndef bivGLWE_CIPHERTEXT_H
#define bivGLWE_CIPHERTEXT_H

#include <stdint.h>

#include "bivariate_polynomial.h"
#include "glwe_key.h"
#include "glwe_params.h"
#include "spqlios_alias.h"

/**
 * @file glwe_ciphertext.h
 *
 * In this header file, we define the structures representing bivariate GLWE ciphertext in both DFT and coefficient forms,
 * encryption operations to them, decryption, allocation and deallocation and accessors to their components.
 */

// bivGLWE PART (begin)

/**
 * @brief A GLWE ciphertext, using bivariate base-2K form
 */
typedef struct glwe_ciphertext
{
	const GLWEParams* params;  ///< GLWE parameters

	/**
	 * @brief Represents a vector of size (k + 1) * l with coefficients that are in \ZnX
	 *
	 * Data layout is limb-major, ie, all the coefficients for a power of \f$ 2^{iK} \f$
	 *  are strored contiguously after the ones from \f$ 2^{(i-1)K} \f$ and before \f$ 2^{(i+1)K} \f$
	 * (see dedicated document for details)
	 *
	 */
	VecBiv* vec;
} GLWECiphertext;

/**
 * @brief Creates a bivGLWE, filled with 0.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return The newly allocated object
 */
GLWECiphertext* new_glwe(const GLWEParams* params_glwe);

/**
 * @brief Deletes a bivGLWE ciphertext, but the bivGLWE parameters.
 *
 * @param glwe The bivGLWE ciphertext.
 */
void delete_glwe(GLWECiphertext* glwe);

/**
 * @brief Encrypts a phase (message + noise) and puts it in result.
 *
 * @param module Additionnal information for backend.
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase message + noise.
 *
 * @retval -1 if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval  0 othwerwise.
 */
int glwe_secret_encrypt_phase(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                              const PolyBiv* phase);

/**
 * @brief Encrypts a univariate message
 *
 * @param module Additionnal information for backend.
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ_rnx  The univariate secret message to encrypt
 *
 * @retval -1 if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval  0 othwerwise.
 */
int glwe_secret_encrypt_rnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                            const PolyUnivRnX* m_univ_rnx);

/**
 * @brief Encrypts a univariate message
 *
 * @param module Additionnal information for backend.
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ_tnx  The univariate secret message to encrypt
 *
 * @retval -1 if an error occurs.
 * @retval  -0 othwerwise.
 */
int glwe_secret_encrypt_tnx(const MODULE* module, GLWECiphertext* result, const GLWESecretKeyDFT* sk_dft,
                            const PolyUnivTnX* m_univ_tnx);

/**
 * @brief Decrypts the ciphertext into the
 * bivariate phase (message + noise) and puts it in result.
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate phase.
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval -1 if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval  0 othwerwise.
 */
int glwe_secret_decrypt(const MODULE* module, PolyBiv* result, const GLWESecretKeyDFT* sk_dft,
                        const GLWECiphertext* glwe);
/**
 * @brief Gives a pointer to the start of a STRIDED polynomial in a GLWECiphertext
 *
 * Critically, it is not a pointer to a contiguous PolyBiv! Due to the memory
 * layout of GLWECiphertext, it is strided, that is, there are other data
 * between the different limbs of the polynomial
 *
 * @param glwe_ct The GLWECiphertext from which to retrieve the start of a bivariate polynomial
 * @param pos The index of the polynomial whose start is to be retrieved
 *
 * @return The start of a strided bivariate polynomial
 *
 */
PolyBiv* glwe_extract_start_poly(const GLWECiphertext* glwe_ct, uint64_t pos);

// bivGLWE IN DFT PART (begin)

/**
 * @brief A GLWE ciphertext in the DFT domain, encoded using bivariate base-2K form
 */
typedef struct glwe_ciphertext_dft
{
	const GLWEParams* params;  ///< GLWE parameters
	VecBivDFT* vec;  ///< Vector in the DFT domain. See GLWECiphertext for layout details, DFT is performed per-limb
} GLWECiphertextDFT;

/**
 * @brief Creates a new empty bivGLWE ciphertext.
 *
 * @param params_glwe The GLWE parameters.
 * @return A pointer to the new object or NULL in case of a failure
 */
GLWECiphertextDFT* new_glwe_dft(const GLWEParams* params_glwe);

/**
 * @brief Deletes a GLWE ciphertext, but not its parameters.
 *
 * @param glwe The GLWE ciphertext.
 */
void delete_glwe_dft(GLWECiphertextDFT* glwe);

/**
 * @brief Gives a pointer to the start of a STRIDED polynomial in a GLWECiphertextDFT
 * Critically, it is not a pointer to a contiguous PolyBivDFT! Due to the memory
 * layout of GLWECiphertextDFT, it is strided, that is, there are other data
 * between the different limbs of the polynomial
 *
 * @param glwe_dft The GLWECiphertext from which to retrieve the start of a bivariate polynomial
 * @param pos The number of the polynomial whose starts is to be retrieved
 *
 * @return PolyBivDFT* The start of a strided bivariate polynomial
 *
 */
PolyBivDFT* glwe_extract_start_poly_dft(const GLWECiphertextDFT* glwe_dft, uint64_t pos);

/**
 * @brief convert a GLWECiphertext into its DFT representation
 *
 * @param module The backend module
 * @param res_dft The resulting glwe ciphertext (in DFT)
 * @param glwe_ct The input coefficient domain GLWE
 *
 */
int glwe_coef_to_dft(const MODULE* module, GLWECiphertextDFT* res_dft, const GLWECiphertext* glwe_ct);

/**
 *
 * @brief Convert a GLWECiphertextDFT into its coefficient (non-DFT) representation
 *
 * @param module The backend module
 * @param res_ct The resulting glwe ciphertext (in coef space)
 * @param glwe_dft The input glwe in the DFT domain
 *
 */
int glwe_dft_to_coef(const MODULE* module, GLWECiphertext* res_ct, const GLWECiphertextDFT* glwe_dft);

#endif  // bivGLWE_CIPHERTEXT_H
