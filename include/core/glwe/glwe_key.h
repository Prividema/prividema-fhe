#ifndef bivGLWE_KEY_H
#define bivGLWE_KEY_H

#include "bivariate_polynomial.h"
#include "glwe_params.h"

// bivGLWE SECRET KEY STRUCTURES

/**
 * @brief Structure that represents a bivGLWE secret key
 *
 */
typedef struct glwe_secret_key
{
	uint64_t nn;      ///< Degree of the polynomials
	uint64_t k;       ///< Number of polynomials of degree nn in key
	VecUniv* values;  ///< Flattened vector of the k polynomials of degree nn
} GLWESecretKey;

/**
 * @brief Structure that represents a bivGLWE secret key in DFT domain
 */
typedef struct glwe_prep_secret_key
{
	uint64_t nn;         ///< Degree of the polynomials
	uint64_t k;          ///< Number of polynomials of degree nn in key
	VecUnivDFT* values;  ///< Flattened vector of the k polynomials of degree nn

} GLWESecretKeyDFT;

// COEF SPACE PART (begin)

/**
 * @brief Creates a GLWE Secret key.
 *
 * @param params_glwe The glwe params for the key
 *
 * @return A pointer to the newly allocated object or NULL if it failed to do so.
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
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int uniform_glwe_secret_key(const MODULE* module, GLWESecretKey* sk, uint64_t nb_bits);

/**
 * @brief Returns a pointer to the k'th polynomial in the secret key
 *
 *
 * @param sk  The secret key
 * @param pos The position to retrieve (from 0 to k-1)
 *
 * @return    A pointer to the beggining of the pos-th polynomial in the key
 */
PolyUniv* glwe_sk_extract_poly(GLWESecretKey* sk, uint64_t pos);

/**
 * @brief Deletes the secret key.
 *
 * @param sk The secret key.
 */
void delete_glwe_secret_key(GLWESecretKey* sk);

// DFT SPACE PART (begin)

/**
 * @brief Creates a GLWE Secret key in the DFT domain.
 *
 * @param params_glwe the glwe params for the key
 *
 * @return A pointer to the newly allocated object or NULL if it failed.
 */
GLWESecretKeyDFT* alloc_glwe_secret_key_dft(GLWEParams* params_glwe);

/**
 * @brief Returns a pointer to the k'th polynomial in the secret key
 *
 *
 * @param sk_dft  The secret key in the DFT domain
 * @param pos The position to retrieve (from 0 to k-1)
 *
 * @return    A pointer to the beggining of the pos-th polynomial in the key
 */
PolyUnivDFT* glwe_sk_extract_poly_dft(const GLWESecretKeyDFT* sk_dft, uint64_t pos);

/**
 * @brief Deletes a secret key in the DFT domain.
 *
 * @param sk_dft The secret key in the DFT domain.
 */
void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft);

////////////////////
// PUBLIC KEY
////////////////////

typedef struct glwe_ciphertext GLWECiphertext;
/**
 * @brief Public bivGLWE key
 *
 * Composed of yy GLWE encryptions of 0
 *
 */

typedef struct glwe_public_key
{
	uint64_t nn;          ///< Degree of the polynomials
	uint64_t k;           ///< Parameter k of the cryptosystem
	int64_t yy;           ///< Number of encryptions of 0 that the public key provides
	GLWECiphertext** pk;  ///< Vector of yy GLWE encryptions of 0
} GLWEPublicKey;

#endif  // GLWE_KEY_H
