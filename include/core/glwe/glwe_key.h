#ifndef bivGLWE_KEY_H
#define bivGLWE_KEY_H

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

//! bivGLWE SECRET KEY STRUCTURES
/**
 * @brief Structure that represents a bivGLWE secret key
 * 
 * @param N The degree of the chosen cyclotomic polynomial.
 */
typedef struct glwe_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUniv** values;
} GLWESecretKey;

typedef struct glwe_prep_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUnivDFT** values; 
} GLWESecretKeyDFT;

//! bivGLWE KEY PART (begin)

/**
 * @brief Creates a new secret key values component.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of Zn[X] polynomials in a secret key.
 * @return PolyUniv**
 */
PolyUniv** new_glwe_secret_key_values(uint64_t N, uint64_t k);

/**
 * @brief Delete a secret key values component.
 *
 * @param values The values of the secret key.
 * @param k The number of Zn[X] polynomials.
 */
void delete_glwe_secret_key_values(PolyUniv** values, uint64_t k);

/**
 * @brief Creates a bivGLWE Secret key.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWESecretKeyDFT*
 */
GLWESecretKey* new_glwe_secret_key(uint64_t N, uint64_t k);

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
 * @brief Delete the secret key.
 *
 * @param sk The secret key.
 */
void delete_glwe_secret_key(GLWESecretKey* sk);


//! bivGLWE IN DFT SPACE PART (begin)

/**
 * @brief Creates a new secret key values component in the DFT domain.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of Zn[X] polynomials in a secret key.
 * @return PolyUnivDFT**
 */
PolyUnivDFT** new_glwe_secret_key_values_dft(uint64_t N, uint64_t k);

/**
 * @brief Delete the values of a secret key in the DFT domain.
 *
 * @param values The values of the secret key in the DFT domain.
 * @param k The number of Zn[X] polynomial in the secret key.
 */
void delete_glwe_secret_key_values_dft(PolyUnivDFT** values, uint64_t k);

/**
 * @brief Creates a bivGLWE Secret key in the DFT domain.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GLWESecretKeyDFT*
 */
GLWESecretKeyDFT* new_glwe_secret_key_dft(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param sk_dft   The secret key in the DFT domain.
 * @param nb_bits  Bit-size of coefficients: coefficients are sampled uniformly
 *                 in the range [-2^nb_bits, 2^nb_bits).
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int uniform_glwe_secret_key_dft(const MODULE* module, GLWESecretKeyDFT* sk_dft, uint64_t nb_bits);

/**
 * @brief Delete the secret key that is in the DFT domain.
 *
 * @param sk_dft The secret key in the DFT domain.
 */
void delete_glwe_secret_key_dft(GLWESecretKeyDFT* sk_dft);

#endif  // bivGLWE_KEY_H