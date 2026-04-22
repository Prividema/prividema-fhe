#ifndef bivGLWE_CIPHERTEXT_H
#define bivGLWE_CIPHERTEXT_H

#include <stdint.h>

#include "bivariate_polynomial.h"
#include "glwe_params.h"
#include "spqlios_alias.h"

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
 * @brief Normalizes a bivGLWE ciphertext.
 *
 * @param module Additionnal information for backend.
 * @param res The result normalized bivGLWE ciphertext.
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int normalize_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* glwe);

/**
 * @brief Adds two bivGLWE ciphertexts.
 *
 * @param res The result bivGLWE ciphertext.
 * @param glwe_lhs The left-hand side bivGLWE ciphertext.
 * @param glwe_rhs The right-hand side bivGLWE ciphertext.
 */
void add_glwe(GLWECiphertext* res, const GLWECiphertext* glwe_lhs, const GLWECiphertext* glwe_rhs);

/**
 * @brief Multiply a bivGLWE ciphertext by a \ZnX polynomial.
 *
 * @param module Additionnal information for backend.
 * @param res The result GLWE ciphertext.
 * @param u The \ZnX polynomial.
 * @param glwe The GLWE ciphertext.
 *
 * @retval -1 if an error occurs
 * @retval 0 otherwise.
 */
int const_mult_glwe(const MODULE* module, GLWECiphertext* res, const PolyUnivDFT* u, const GLWECiphertext* glwe);

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
 * @brief The number of coefficient in a bivariate GLWE ciphertext in the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return The number of coefficient in a bivariate GLWE ciphertext in the DFT domain.
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * \ZnX, due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t glwe_coef_number_dft(const GLWEParams* params_glwe);

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
 * @brief Adds two bivGLWE ciphertexts.
 *
 * @param res_dft The result bivGLWE ciphertext in the DFT domain.
 * @param glwe_lhs_dft The left-hand side bivGLWE ciphertext in the DFT domain.
 * @param glwe_rhs_dft The right-hand side bivGLWE ciphertext in the DFT domain.
 */
void add_glwe_dft(GLWECiphertextDFT* res_dft, const GLWECiphertextDFT* glwe_lhs_dft,
                  const GLWECiphertextDFT* glwe_rhs_dft);

/**
 * @brief Multiply a bivGLWE ciphertext by a \ZnX polynomial in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param res_dft The result bivGLWE ciphertext in the DFT domain.
 * @param u The \ZnX polynomial.
 * @param glwe_dft The bivGLWE ciphertext in the DFT domain.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int const_mult_glwe_dft(const MODULE* module, GLWECiphertextDFT* res_dft, const PolyUnivDFT* u,
                        const GLWECiphertextDFT* glwe_dft);

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

// COMMON PART (begin)

/**
 * @brief Compute the polynomial product of c and d, component-wise in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param res_dft The result in the DFT domain.
 * @param res_size The result's size.
 * @param c_dft The left-hand side polynomial in the DFT domain .
 * @param c_size The left-hand size of c_dft.
 * @param d_dft The right-hand side polynomial in the DFT domain.
 * @param d_size The right-hand size of c_dft.
 *
 * @remark `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the
 * size of res_dft.
 */
void mult_vec_znx_dft(const MODULE* module, double* res_dft, int64_t res_size, const double* c_dft, int64_t c_size,
                      const double* d_dft, int64_t d_size);

#endif  // bivGLWE_CIPHERTEXT_H
