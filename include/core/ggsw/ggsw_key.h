#ifndef GGSW_KEY_H
#define GGSW_KEY_H

/**
 * @file ggsw_key.h
 *
 * In this header file, we define the structure representing GGSW secret key in both DFT and iDFT forms.
 */

#include "ggsw_ciphertext.h"
#include "glwe_key.h"

typedef struct ggsw_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUniv** values;  // The key itself. Vector of size k.
} GGSWSecretKey;

typedef struct ggsw_prep_secret_key
{
	uint64_t N;
	uint64_t k;
	PolyUnivDFT** values;  // The key itself in the DFT domain. Vector of size k, each element is in the DFT domain
	void* data;			   // TODO : What does data contain if values exists ????
} GGSWSecretKeyDFT;

// =============================================
// |                                           |
// |      			 iDFT Form       		   |
// |                                           |
// =============================================

/**
 * @brief Creates a new secret key values component.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of \f$\mathbb{Z}_n[X]\f$ polynomials in a secret key.
 * 
 * @retval • `NULL` if malloc failed inside the function.
 * @retval • Allocated secret key's values. 
 */
PolyUniv** new_ggsw_secret_key_values(uint64_t N, uint64_t k);

/**
 * @brief Deletes a secret key values component.
 *
 * @param values The values of the secret key.
 * @param k 	 The number of \f$\mathbb{Z}_n[X]\f$ polynomials.
 */
void delete_ggsw_secret_key_values(PolyUniv** values, uint64_t k);

/**
 * @brief Creates a GGSW secret key.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of \f$\mathbb{Z}_n[X]\f$ polynomials in the secret key.
 * 
 * @retval • `NULL` if malloc failed inside the function.
 * @retval • Allocated secret key, filled with 0.
 */
GGSWSecretKey* new_ggsw_secret_key(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly.
 *
 * @param module  Additionnal information for backend.
 * @param sk 	  The secret key.
 * @param nb_bits The exponent of the range = [-2^nb_bits, 2^nb_bits).
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 */
int uniform_ggsw_secret_key(const MODULE* module, GGSWSecretKey* sk, uint64_t nb_bits);

/**
 * @brief Deletes the secret key.
 *
 * @param sk The secret key.
 */
void delete_ggsw_secret_key(GGSWSecretKey* sk);


// =============================================
// |                                           |
// |      			 DFT Form       		   |
// |                                           |
// =============================================


/**
 * @brief Creates a new secret key values component in the DFT domain.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of \f$\mathbb{Z}_n[X]\f$ polynomials in a secret key.
 * 
 * @retval • `NULL` if malloc failed inside the function.
 * @retval • Allocated secret key's values.
 */
PolyUnivDFT** new_ggsw_secret_key_values_dft(uint64_t N, uint64_t k);

/**
 * @brief Deletes the values of a secret key in the DFT domain.
 *
 * @param values The values of the secret key in the DFT domain.
 * @param k 	 The number of \f$\mathbb{Z}_n[X]\f$ polynomials in a secret key.
 */
void delete_ggsw_secret_key_values_dft(PolyUnivDFT** values, uint64_t k);

/**
 * @brief Creates a GGSW secret key in the DFT domain.
 *
 * @param N The degree of the chosen cyclotomic polynomial.
 * @param k The number of \f$\mathbb{Z}_n[X]\f$ polynomials in the secret key.
 * 
 * @return A Pointer to the GGSW secret key in the DFT domain.
 */
GGSWSecretKeyDFT* new_ggsw_secret_key_dft(uint64_t N, uint64_t k);

/**
 * @brief Draws a secret key uniformly in the DFT domain.
 *
 * @param module   Additionnal information for backend.
 * @param sk_dft   The secret key in the DFT domain.
 * @param nb_bits  Bit-size of coefficients: coefficients are sampled uniformly
 *                 in the range [-2^nb_bits, 2^nb_bits).
 *
 * @retval • `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval • `0` otherwise.
 */
int uniform_ggsw_secret_key_dft(const MODULE* module, GGSWSecretKeyDFT* sk_dft, uint64_t nb_bits);

/**
 * @brief Deletes the secret key that is in the DFT domain.
 *
 * @param sk_dft A Pointer to the secret key in the DFT domain.
 */
void delete_ggsw_secret_key_dft(GGSWSecretKeyDFT* sk_dft);

#endif  // GGSW_KEY_H
