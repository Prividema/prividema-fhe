#ifndef bivGLWE_CIPHERTEXT_H
#define bivGLWE_CIPHERTEXT_H

#include "bivariate_polynomial.h"

//! bivGLWE PART (begin)

typedef struct glwe_ciphertext
{
	const GLWECtParams* params;  // bivGLWE parameters
	VecBiv* vec;                 // Represents a vector of size (k + 1) * l with coefficients that are in Zn[X]
} GLWECiphertext;

/**
 * @brief Return the number of coefficient in a bivariate bivGLWE ciphertext.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return int64_t
 */
uint64_t glwe_coef_number(const GLWECtParams* params_glwe);

/**
 * @brief Creates a bivGLWE, filled with 0.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return GLWECiphertext*
 */
GLWECiphertext* new_glwe(const GLWECtParams* params_glwe);

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
 * @param result The result normalized bivGLWE ciphertext.
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int normalize_glwe(const MODULE* module, GLWECiphertext* res, const GLWECiphertext* glwe);

/**
 * @brief Adds two bivGLWE ciphertexts.
 *
 * @param result The result bivGLWE ciphertext.
 * @param glwe_lhs The left-hand side bivGLWE ciphertext.
 * @param glwe_rhs The right-hand side bivGLWE ciphertext.
 */
void add_glwe(GLWECiphertext* res, const GLWECiphertext* glwe_lhs, const GLWECiphertext* glwe_rhs);

/**
 * @brief Multiply a bivGLWE ciphertext by a Zn[X] polynomial.
 *
 * @param module Additionnal information for backend.
 * @param result The result bivGLWE ciphertext.
 * @param u The Zn[X] polynomial.
 * @param glwe The bivGLWE ciphertext.
 * @param do_normalization The function normalizes the bivGLWE ciphertext if and only if do_normalization = 1.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int const_mult_glwe(const MODULE* module, GLWECiphertext* res, const PolyUnivDFT* u, const GLWECiphertext* glwe,
                    int do_normalization);

//! bivGLWE IN DFT PART (begin)

typedef struct glwe_ciphertext_dft
{
	const GLWECtParams* params;  // bivGLWE parameters
	VecBivDFT* vec;              // Vector in the DFT
} GLWECiphertextDFT;

/**
 * @brief The number of coefficient in a bivariate bivGLWE ciphertext in the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return int64_t
 *
 * @note The number of independent coefficients of a polynomial in the DFT domain is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t glwe_coef_number_dft(const GLWECtParams* params_glwe);

/**
 * @brief Creates a new empty bivGLWE ciphertext.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return GLWECiphertextDFT*
 */
GLWECiphertextDFT* new_glwe_dft(const GLWECtParams* params_glwe);

/**
 * @brief Deletes a bivGLWE ciphertext, but not the parameters.
 *
 * @param glwe The bivGLWE ciphertext.
 */
void delete_glwe_dft(GLWECiphertextDFT* glwe);

// TODO test normalize_glwe_dft
/**
 * @brief Normalizes a bivGLWE ciphertext in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result normalized bivGLWE ciphertext in the DFT domain.
 * @param glwe_dft The bivGLWE ciphertext in the DFT domain.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int normalize_glwe_dft(const MODULE* module, GLWECiphertextDFT* result_dft, const GLWECiphertextDFT* glwe_dft);

/**
 * @brief Adds two bivGLWE ciphertext.
 *
 * @param res_dft The result bivGLWE ciphertext in the DFT domain.
 * @param glwe_lhs_dft The left-hand side bivGLWE ciphertext in the DFT domain.
 * @param glwe_rhs_dft The right-hand side bivGLWE ciphertext in the DFT domain.
 */
void add_glwe_dft(GLWECiphertextDFT* res_dft, const GLWECiphertextDFT* glwe_lhs_dft,
                  const GLWECiphertextDFT* glwe_rhs_dft);

/**
 * @brief Multiply a bivGLWE ciphertext by a Zn[X] polynomial in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param res_dft The result bivGLWE ciphertext in the DFT domain.
 * @param u The Zn[X] polynomial.
 * @param glwe_dft The bivGLWE ciphertext in the DFT domain.
 * @param do_normalization the function normalizes the bivGLWE ciphertext if and only if do_normalization = 1.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int const_mult_glwe_dft(const MODULE* module, GLWECiphertextDFT* res_dft, const PolyUnivDFT* u,
                        const GLWECiphertextDFT* glwe_dft, int do_normalization);

//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGLWE ciphertext, in the DFT domain & out of the DFT domain.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return int64_t
 *
 * @note The size of a bivGLWE ciphertext is the same in and out of the DFT domain.
 */
uint64_t glwe_size(const GLWECtParams* params_glwe);

/**
 * @brief The number of bytes needed to store a bivGLWE ciphertext.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return int64_t
 *
 * @note The number of bytes needed to store a bivGLWE ciphertext, is the same in and out of the DFT domain.
 */
uint64_t glwe_bytes(const GLWECtParams* params_glwe);

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
 * @note `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the
 * size of res_dft.
 */
void mult_vec_znx_dft(const MODULE* module, double* res_dft, int64_t res_size, const double* c_dft, int64_t c_size,
                      const double* d_dft, int64_t d_size);

#endif  // bivGLWE_CIPHERTEXT_H
