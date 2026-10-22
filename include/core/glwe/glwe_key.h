#ifndef GLWE_KEY_H
#define GLWE_KEY_H

#include "backend.h"
#include "backend_types.h"
#include "glwe_params.h"
#include "univariate_polynomial.h"

/**
 *
 * @file glwe_key.h
 *
 * @brief Definition and operations for GLWE/GLWEGadget/GGSW secret keys.
 */

// =============================================
// |                                           |
// |                Secret Key                 |
// |                                           |
// =============================================

/**
 * @brief Structure that represents a GLWE secret key,
 * which is also the same as a GGSW or GLWEGadget secret key
 */
typedef struct glwe_secret_key
{
	uint64_t nn;      ///< Ring dimension (also the number of coefficients in one polynomial)
	uint64_t k;       ///< Number of polynomials in the key
	int64_t* values;  ///< Flattened vector containing the data
} GLWESecretKey;

/**
 * @brief Creates a GLWE Secret key. This function allocates the memory.
 *
 * @param params_glwe The GLWE params for the key.
 *
 * @return A pointer to the newly allocated object or NULL if it failed to do so.
 */
GLWESecretKey* new_glwe_secret_key(GLWEParams* params_glwe);

/**
 * @brief Deletes the secret key. This function frees the memory.
 *
 * @param sk The secret key to delete.
 */
void delete_glwe_secret_key(GLWESecretKey* sk);

/**
 * @brief Draws a secret key uniformly.
 *
 * @param module  The backend instance.
 * @param sk 	  The secret key.
 * @param nb_bits Bit-size of coefficients: coefficients are sampled uniformly
 *                in the range [-2^(nb_bits-1), 2^(nb_bits-1)).
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int glwe_sk_random_uniform(const PvdaBackend* module, GLWESecretKey* sk, uint64_t nb_bits);

/**
 * @brief Draws a binary secret key.
 *
 * @param module The backend instance.
 * @param sk 	 The secret key.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int glwe_sk_random_binary(const PvdaBackend* module, GLWESecretKey* sk);

/**
 * @brief Draws a ternary secret key.
 *
 * @param module The backend instance.
 * @param sk 	 The secret key.
 *
 * @retval -1 if an error occurs.
 * @retval 0 otherwise.
 */
int glwe_sk_random_ternary(const PvdaBackend* module, GLWESecretKey* sk);

/**
 * @brief Returns a view of the pos'th polynomial in the secret key.
 *
 * @param sk  The secret key
 * @param pos The position to retrieve (from 0 to k-1)
 *
 * @return A pointer to the N coefficients of the pos'th polynomial, inside sk.
 */
PolyUniv* glwe_sk_extract_poly(GLWESecretKey* sk, uint64_t pos);

// =============================================
// |                                           |
// |           Secret Key Prepared             |
// |                                           |
// =============================================

/**
 * @brief Structure that represents a prepared GLWE secret key.
 * 
 * We differentiate "prepared" and "in DFT" as in spqlios
 * "DFT" is a general-purpose representation while
 * "prepared" is an operand reshaped for one operation (SVP or VMP).
 * 
 * Here the key is prepared for SVP.
 */
typedef struct glwe_secret_key_prepared
{
	const PvdaBackend* module; 		///< The backend instance
	uint64_t nn;           	   		///< Ring dimension
	uint64_t k;            	   		///< Number of polynomials in the key
	PolyUniv* values_coef;         	///< The k key polynomials in coefficient space, contiguous (k * nn coefficients)
	BE_PolySvpPrep** values_prep;  	///< The k key polynomials as SVP operands (k separate backend objects)
} GLWESecretKeyPrepared;

/**
 * @brief Creates a prepared GLWE secret key. This function allocates the memory.
 *
 * @param params_glwe the GLWE params for the key.
 *
 * @return A pointer to the newly allocated object or NULL if it failed.
 */
GLWESecretKeyPrepared* new_glwe_secret_key_prepared(const PvdaBackend* module, GLWEParams* params_glwe);

/**
 * @brief Deletes a prepared secret key
 *
 * @param sk_prep The prepared secret key
 */
void delete_glwe_secret_key_prepared(GLWESecretKeyPrepared* sk_prep);

/**
 * @brief Returns the pos'th component of the prepared secret key
 *
 * @param sk_prep The prepared secret key
 * @param pos     The position to retrieve (from 0 to k-1)
 *
 * @return The pos'th key polynomial, prepared as the scalar operand of scalar-vector products.
 *         It is owned by sk_prep: do not delete it.
 */
BE_PolySvpPrep* glwe_prepared_sk_extract_poly_prep(const GLWESecretKeyPrepared* sk_prep, uint64_t pos);

/**
 * @brief Returns a view of the pos'th polynomial in the prepared secret key.
 *
 *
 * @param sk_prep  The prepared secret key.
 * @param pos The position to retrieve (from 0 to k-1).
 *
 * @return A pointer to the beginning of the pos-th polynomial in the key, in coefficient form.
 */
PolyUniv* glwe_sk_prepared_extract_poly_coefs(const GLWESecretKeyPrepared* sk_prep, uint64_t pos);


// =============================================
// |                                           |
// |           Secret Key Transforms           |
// |                                           |
// =============================================

/**
 * @brief Prepares a secret key: copies its coefficients and prepares each component
 * as the scalar operand of scalar-vector products
 *
 * @param module      The backend instance (must be the one sk_prep was allocated with)
 * @param sk_prep     The resulting prepared GLWE secret key
 * @param sk          The secret key in coefficient space
 *
 * @retval -1 if an error occurs
 *
 * @retval 0 otherwise
 */
int glwe_sk_prepare(const PvdaBackend* module, GLWESecretKeyPrepared* sk_prep, const GLWESecretKey* sk);

// =============================================
// |                                           |
// |                Public Key                 |
// |                                           |
// =============================================

// We have to declare GLWE ciphertext here as a public key is just GLWE encryptions of 0.
// But glwe_ciphertext.h imports this file, so to avoid include loops the structure is declared here.
typedef struct glwe_ciphertext GLWECiphertext;

/**
 * @brief Public bivariate GLWE key
 *
 * Composed of yy bivariate GLWE encryptions of 0
 *
 */
typedef struct glwe_public_key
{
	uint64_t nn;          ///< Ring dimension
	uint64_t k;           ///< Number of polynomials in the key
	int64_t size_pk;      ///< Number of encryptions of 0 that the public key provides
	GLWECiphertext** pk;  ///< Vector of GLWE encryptions of 0
} GLWEPublicKey;

#endif  // GLWE_KEY_H
