#ifndef bivGLWE_KEY_H
#define bivGLWE_KEY_H

#include "bivariate_polynomial.h"
#include "glwe_ciphertext.h"
#include "glwe_params.h"

//! bivGLWE SECRET KEY STRUCTURES
/**
 * @brief Structure that represents a bivGLWE secret key
 *
 * @param nn The degree of the chosen cyclotomic polynomial.
 */
typedef struct glwe_secret_key
{
	uint64_t nn;
	uint64_t k;
	PolyUniv* values;
} GLWESecretKey;

typedef struct glwe_prep_secret_key
{
	uint64_t nn;
	uint64_t k;
	PolyUnivDFT* values;
} GLWESecretKeyDFT;

//! bivGLWE KEY PART (begin)

/**
 * @brief Creates a GLWE Secret key.
 *
 * @params the glwe params for the key
 *
 * @return GLWESecretKeyDFT*
 */
GLWESecretKey* alloc_glwe_secret_key(GLWEParams* params_glwe);

/**
 * @brief Draws a secret key uniformly.
 *
 * @param module Additionnal information for backend.
 * @param sk The secret key.
 * @param nb_bits  Bit-size of coefficients: coefficients are sampled uniformly
 *                 in the range [-2^nb_bits, 2^nb_bits).
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int uniform_glwe_secret_key(const MODULE* module, GLWESecretKey* sk, uint64_t nb_bits);

/**
 * Returns a pointer to the k'th polynomial in the secret key
 *
 *
 * @param sk  The secret key
 * @param pos The position to retrieve (from 0 to k-1)
 *
 * @return    A pointer to the beggining of the pos-th polynomial in the key
 */
PolyUniv* glwe_sk_extract_poly(GLWESecretKey* sk, uint64_t pos);

/**
 * @brief Delete the secret key.
 *
 * @param sk The secret key.
 */
void delete_glwe_secret_key(GLWESecretKey* sk);

//! bivGLWE IN DFT SPACE PART (begin)

/**
 * @brief Creates a GLWE Secret key in the DFT domain.
 *
 * @param params_glwe the glwe params for the key
 *
 * @return GLWESecretKeyDFT*
 */
GLWESecretKeyDFT* alloc_glwe_secret_key_dft(GLWEParams* params_glwe);

/**
 * Returns a pointer to the k'th polynomial in the secret key
 *
 *
 * @param sk  The secret key
 * @param pos The position to retrieve (from 0 to k-1)
 *
 * @return    A pointer to the beggining of the pos-th polynomial in the key
 */
PolyUnivDFT* glwe_sk_extract_poly_dft(const GLWESecretKeyDFT* sk_dft, uint64_t pos);

/**
 * @brief Delete the secret key that is in the DFT domain.
 *
 * @param sk_dft The secret key in the DFT domain.
 */
void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft);

typedef struct glwe_public_key
{
	uint64_t nn;
	uint64_t k;
	int64_t yy;
	GLWECiphertext** pk;  // vector of Y element (A1,...,Ak, B)
} GLWEPublicKey;

#endif  // GLWE_KEY_H
