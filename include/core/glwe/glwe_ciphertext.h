#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>

#include "bivariate_polynomial.h"
#include "glwe_ct_params.h"
#include "spqlios_alias.h"
#include "vec_znx_arithmetic_private.h"

//! GLWE PART (begin)

typedef struct glwe_ciphertext
{
	GLWECtParams* params;  // GLWE parameters
	VecBiv* vec;           // Represents a vector of size (k + 1) * l with coefficients that are in Zn[X]
} GLWECiphertext;

/**
 * @brief Return the number of coefficient in a bivariate GLWE ciphertext.
 *
 * @param params The GLWE parameters.
 * @return int64_t
 */
uint64_t glwe_coef_number(GLWECtParams* params);

/**
 * @brief Creates a bivGLWE, filled with 0.
 *
 * @param params The GLWE parameters.
 * @return GLWECiphertext*
 */
GLWECiphertext* new_glwe(GLWECtParams* params);

/**
 * @brief Deletes a GLWE ciphertext, but the GLWE parameters.
 *
 * @param ct The GLWE ciphertext.
 */
void delete_glwe(GLWECiphertext* ct);

/**
 * @brief Normalizes a GLWE ciphertext.
 *
 * @param module
 * @param res The result normalized GLWE ciphertext.
 * @param ct_glwe The GLWE ciphertext.
 */
void normalize_glwe(MODULE* module, GLWECiphertext* res, GLWECiphertext* ct_glwe);

/**
 * @brief Adds two GLWE ciphertexts.
 *
 * @param res The result GLWE ciphertext.
 * @param ct1 The left-hand side GLWE ciphertext.
 * @param ct2 The right-hand side GLWE ciphertext.
 */
void add_glwe(GLWECiphertext* res, GLWECiphertext* ct1, GLWECiphertext* ct2);

/**
 * @brief Multiply a GLWE ciphertext by a Zn[X] polynomial.
 *
 * @param module
 * @param res The result GLWE ciphertext.
 * @param u The Zn[X] polynomial.
 * @param ct The GLWE ciphertext.
 * @param do_normalization The function normalizes the GLWE ciphertext if and only if do_normalization = 1.
 */
int const_mult_glwe(MODULE* module, GLWECiphertext* res, PolyUnivDFT* u, GLWECiphertext* ct, int do_normalization);

//! GLWE IN DFT PART (begin)

typedef struct glwe_ciphertext_dft
{
	GLWECtParams* params;  // GLWE parameters
	VecBivDFT* vec;       // Prepared vector
} GLWECiphertextDFT;

/**
 * @brief The number of coefficient in a bivariate GLWE ciphertext in DFT space.
 *
 * @param params The GLWE parameters.
 * @return int64_t
 *
 * @note The number of independent coefficients of a polynomial in DFT space is half the number of coefficients in
 * Zn[X], due to conjugate symmetry when the polynomial has real (or integer) coefficients.
 */
uint64_t glwe_coef_number_dft(GLWECtParams* params);

/**
 * @brief Creates a new empty GLWE ciphertext.
 *
 * @param params The GLWE parameters.
 * @return GLWECiphertextDFT*
 */
GLWECiphertextDFT* new_glwe_dft(GLWECtParams* params);

/**
 * @brief Deletes a GLWE ciphertext, but not the parameters.
 *
 * @param ct The GLWE ciphertext.
 */
void delete_glwe_dft(GLWECiphertextDFT* ct);

/**
 * @brief Adds two GLWE ciphertext.
 *
 * @param res_dft The result GLWE ciphertext in DFT space.
 * @param ct1_dft The left-hand side GLWE ciphertext in DFT space.
 * @param ct2_dft The right-hand side GLWE ciphertext in DFT space.
 */
void add_glwe_dft(GLWECiphertextDFT* res_dft, GLWECiphertextDFT* ct1_dft, GLWECiphertextDFT* ct2_dft);

/**
 * @brief Multiply a GLWE ciphertext by a Zn[X] polynomial in DFT space.
 *
 * @param module
 * @param res_dft The result GLWE ciphertext in DFT space.
 * @param u The Zn[X] polynomial.
 * @param ct_dft The GLWE ciphertext in DFT space.
 * @param do_normalization the function normalizes the GLWE ciphertext if and only if do_normalization = 1.
 */
int const_mult_glwe_dft(MODULE* module, GLWECiphertextDFT* res_dft, PolyUnivDFT* u, GLWECiphertextDFT* ct_dft,
                         int do_normalization);

//! COMMON PART (begin)

/**
 * @brief Return the size of a bivGLWE ciphertext, in DFT space & out of DFT space.
 *
 * @param params The GLWE parameters.
 * @return int64_t
 *
 * @note The size of a bivGLWE ciphertext is the same in and out of DFT space.
 */
uint64_t glwe_size(GLWECtParams* params);

/**
 * @brief The number of bytes needed to store a bivGLWE ciphertext.
 *
 * @param params The GLWE parameters.
 * @return int64_t
 *
 * @note The number of bytes needed to store a bivGLWE ciphertext, is the same in and out of DFT space.
 */
uint64_t glwe_bytes(GLWECtParams* params);

/**
 * @brief Compute the polynomial product of c and d, component-wise in DFT space.
 *
 * @param module The module stocking the degree N.
 * @param res_dft The result in DFT space.
 * @param res_size The result's size.
 * @param c_dft The left-hand side polynomial in DFT space .
 * @param c_size The left-hand size of c_dft.
 * @param d_dft The right-hand side polynomial in DFT space.
 * @param d_size The right-hand size of c_dft.
 *
 * @note `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the
 * size of res_dft.
 */
void mult_vec_znx_dft(const MODULE* module, double* res_dft, int64_t res_size, double* c_dft, int64_t c_size,
                      double* d_dft, int64_t d_size);

#endif  // GLWE_CIPHERTEXT_H