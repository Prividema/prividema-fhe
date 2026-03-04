#ifndef GLWE_KEY_H
#define GLWE_KEY_H

#include "glwe_ciphertext.h"

// TODO : To be Defined (spqlios).
/**
 * @brief obtain a info for:
 *  - the dimension N (or the complex dimension m=N/2)
 *  - any moduleuted fft or ntt items
 *  - the hardware (avx, arm64, x86, ...)
 */
typedef struct core
{
} Core;

//! GLWE SECRET KEY STRUCTURES
/**
 * @brief Structure that represents a GLWE secret key
 * 
 * @param N Le degré du 2N-ème polynome cyclotomic
 */
typedef struct glwe_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUniv** values;
	void* data;
} GLWESecretKey;

typedef struct glwe_prep_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUnivDFT** values; 
	void* data;
} GLWESecretKeyDFT;

//! GLWE KEY PART (begin)

/**
 * @brief Creates a new secret key values component.
 *
 * @param N The dimension of the module Zn[X].
 * @param k The number of Zn[X] polynomials in a secret key.
 * @return PolyUniv**
 */
PolyUniv** new_glwe_secret_key_values(uint64_t N, uint64_t k);

/**
 * @brief Creates a new secret key values component, where each Zn[X] poly is uniformly drawn.
 *
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param N The dimension of the module Zn[X].
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * @return PolyUniv**
 */
PolyUniv** new_uniform_glwe_secret_key_values(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete a secret key values component.
 *
 * @param values The values of the secret key.
 * @param k The number of Zn[X] polynomials.
 */
void delete_glwe_secret_key_values(PolyUniv** values, uint64_t k);

/**
 * @brief Creates a GLWE Secret key.
 *
 * @param values The values of the secret key.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWESecretKeyDFT*
 */
GLWESecretKey* new_glwe_secret_key(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly.
 *
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * @return GLWESecretKeyDFT*
 */
GLWESecretKey* new_uniform_glwe_secret_key(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the secret key.
 *
 * @param sk The secret key.
 */
void delete_glwe_secret_key(GLWESecretKey* sk);

/**
 * @brief Computes the secret key out of DFT space.
 *
 * @param sk_dft The secret key in the DFT domain.
 * @return GLWESecretKey*
 */
GLWESecretKey* transform_glwe_secret_key_dft_to_not_dft(const MODULE* module, const GLWESecretKeyDFT* sk_dft);

/**
 * @brief Computes the values of the secret key out of DFT space.
 *
 * @param values_dft The values of the secret key in the DFT domain.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUniv**
 */
PolyUniv** transform_glwe_secret_key_values_dft_to_not_dft(const MODULE* module, const PolyUnivDFT** values_dft, uint64_t k);

//! GLWE IN DFT SPACE PART (begin)

/**
 * @brief Creates a new secret key values component in the DFT domain.
 *
 * @param N The dimension of the module Zn[X].
 * @param k The number of Zn[X] polynomials in a secret key.
 * @return PolyUnivDFT**
 */
PolyUnivDFT** new_glwe_secret_key_values_dft(uint64_t N, uint64_t k);

/**
 * @brief Creates a new secret key values component in the DFT domain, where each Zn[X] poly is uniformly drawn.
 *
 * @param k The number of Zn[X] polynomials in a secret key.
 * @param N The dimension of the module Zn[X].
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 * @return PolyUnivDFT**
 */
PolyUnivDFT** new_uniform_glwe_secret_key_values_dft(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the values of a secret key in the DFT domain.
 *
 * @param values The values of the secret key in the DFT domain.
 */
void delete_glwe_secret_key_values_dft(PolyUnivDFT** values, uint64_t k);

/**
 * @brief Creates a GLWE Secret key in the DFT domain.
 *
 * @param values The values of the secret key in the DFT domain.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWESecretKeyDFT*
 */
GLWESecretKeyDFT* new_glwe_secret_key_dft(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly in the DFT domain.
 *
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWESecretKeyDFT*
 */
GLWESecretKeyDFT* new_uniform_glwe_secret_key_dft(const MODULE* module, uint64_t k, uint64_t nb_bits);

/**
 * @brief Delete the secret key that is in the DFT domain.
 *
 * @param sk_dft The secret key in the DFT domain.
 */
void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft);

/**
 * @brief Computes the secret key in the DFT domain.
 *
 * @param sk The secret key out of DFT space.
 * @return GLWESecretKey*
 */
GLWESecretKeyDFT* transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, const GLWESecretKey* sk);

/**
 * @brief Computes the values of the secret key in the DFT domain.
 *
 * @param values The values of the secret key out of DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUnivDFT**
 */
PolyUnivDFT** transform_glwe_secret_key_values_not_dft_to_dft(const MODULE* module, const PolyUniv** values, uint64_t k);

typedef struct glwe_public_key
{
	uint64_t N;
	uint64_t k;
	int64_t Y;
	GLWECiphertext** pk;  // vector of Y element (A1,...,Ak, B)
	void* data;
} GLWEPublicKey;

// generation secret key
void glwe_secret_key_gen(const Core* core, uint64_t lambda, GLWESecretKey* s);

// generation public key
void glwe_public_key_gen(const Core* core, uint64_t lambda, GLWESecretKey* s, GLWEPublicKey* pk);

#endif  // GLWE_KEY_H