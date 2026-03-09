#ifndef GGSW_KEY_H
#define GGSW_KEY_H

#include "ggsw_ciphertext.h"
#include "glwe_key.h"

typedef struct ggsw_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUniv** values;  // The key itself.
} GGSWSecretKey;

typedef struct ggsw_prep_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUnivDFT** values;  // vec of size k, each element is prepared vec
	void* data;
} GGSWSecretKeyDFT;

//! GGSW KEY PART (begin)

/**
 * @brief Creates a new secret key values component.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of Zn[X] polynomials in a secret key.
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key's values. 
 */
PolyUniv** new_ggsw_secret_key_values(uint64_t N, uint64_t k);

/**
 * @brief Creates a new secret key values component, where each Zn[X] poly is uniformly drawn.
 *
 * @param module Additionnal information for backend.
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key's uniformly drawn values. 
 */
PolyUniv** new_uniform_ggsw_secret_key_values(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete a secret key values component.
 *
 * @param values The values of the secret key.
 * @param k The number of Zn[X] polynomials.
 */
void delete_ggsw_secret_key_values(PolyUniv** values, uint64_t k);

/**
 * @brief Creates a GGSW secret key.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of Zn[X] polynomial in the secret key.
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key, filled with 0.
 */
GGSWSecretKey* new_ggsw_secret_key(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly.
 *
 * @param module Additionnal information for backend.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key, filled with uniformly drawn value.
 */
GGSWSecretKey* new_uniform_ggsw_secret_key(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the secret key.
 *
 * @param sk The secret key.
 */
void delete_ggsw_secret_key(GGSWSecretKey* sk);

/**
 * @brief Computes the input secret key out of the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param sk_dft The secret key in the DFT domain.
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key. It is the input secret key int the DFT domain, computed of the DFT domain.
 */
GGSWSecretKey* transform_ggsw_secret_key_dft_to_not_dft(const MODULE* module, const GGSWSecretKeyDFT* sk_dft);

/**
 * @brief Computes the values of the input secret key out of the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param values_dft The values of the secret key in the DFT domain.
 * @param k The number of Zn[X] polynomial in the secret key.
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key's values. It is the input secret key's values in the DFT domain, computed out of the DFT domain.
 */
PolyUniv** transform_ggsw_secret_key_values_dft_to_not_dft(const MODULE* module, const PolyUnivDFT** values_dft, uint64_t k);

/**
 * @brief Transforms a GGSW secret key to a GLWE secret key.
 *
 * @param sk_ggsw The GGSW secret key
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated GLWE secret key. The GLWE sk's values equal the GGSW sk's values.
 */
GLWESecretKey* transform_ggsw_secret_key_to_glwe_secret_key(const GGSWSecretKey* sk_ggsw);

//! GGSW IN DFT SPACE PART (begin)

/**
 * @brief Creates a new secret key values component in the DFT domain.
 *
 * @param N The dimension of the module Zn[X].
 * @param k The number of Zn[X] polynomials in a secret key.
 * 
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key's values.
 */
PolyUnivDFT** new_ggsw_secret_key_values_dft(uint64_t N, uint64_t k);

/**
 * @brief Creates a new secret key values component in the DFT domain, where each Zn[X] poly is uniformly drawn.
 *
 * @param module Additionnal information for backend.
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval `NULL` if malloc failed inside the function.
 * @retval Allocated secret key's values.
 */
PolyUnivDFT** new_uniform_ggsw_secret_key_values_dft(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the values of a secret key in the DFT domain.
 *
 * @param values The values of the secret key in the DFT domain.
 * @param k The number of Zn[X] polynomials in a secret key.
 */
void delete_ggsw_secret_key_values_dft(PolyUnivDFT** values, uint64_t k);

/**
 * @brief Creates a GGSW secret key in the DFT domain.
 *
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GGSWSecretKeyDFT*
 */
GGSWSecretKeyDFT* new_ggsw_secret_key_dft(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly in the DFT domain.
 *
 * @param module   TODO
 * @param k        Number of polynomials in the secret key (dimension of the key).
 * @param nb_bits  Bit-size of coefficients: coefficients are sampled uniformly
 *                 in the range [-2^nb_bits, 2^nb_bits).
 *
 * @return Pointer to a newly allocated GGSWSecretKeyDFT, or NULL on failure.
 */
GGSWSecretKeyDFT* new_uniform_ggsw_secret_key_dft(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the secret key that is in the DFT domain.
 *
 * @param sk_dft The secret key in the DFT domain.
 */
void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft);

/**
 * @brief Computes the secret key in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param sk The secret key out of the DFT domain.
 * @return GGSWSecretKey*
 */
GGSWSecretKeyDFT* transform_ggsw_secret_key_not_dft_to_dft(const MODULE* module, const GGSWSecretKey* sk);

/**
 * @brief Computes the values of the secret key in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param values The values of the secret key out of the DFT domain.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUnivDFT**
 */
PolyUnivDFT** transform_ggsw_secret_key_values_not_dft_to_dft(const MODULE* module, const PolyUniv** values, uint64_t k);

/**
 * @brief Transforms a GGSW secret key in the DFT domain to a GLWE secret key in the DFT domain.
 *
 * @param sk_ggsw_dft The GGSW secret key in the DFT domain.
 * @return GLWESecretKey*
 */
GLWESecretKeyDFT* transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(const GGSWSecretKeyDFT* sk_ggsw_dft);

typedef struct ggsw_public_key
{
	uint32_t size;       // number of ciphertexts in public key
	GGSWCiphertext* pk;  // Public key is multiple encryptions of 0
} GGSWPublicKey;

#endif  // GGSW_KEY_H
