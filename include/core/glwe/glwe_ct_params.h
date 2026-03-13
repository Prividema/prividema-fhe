#ifndef bivGLWE_CT_PARAMS_H
#define bivGLWE_CT_PARAMS_H

#include <stdint.h>
#include <stdlib.h>

typedef struct glwe_ct_params
{
	uint64_t N;      // The cyclotomic polynomial degree, in the code m = N / 2, ie the dimension of the DFT domain
	uint64_t k;      // Number of secret key terms, k=1 for RLWE
	uint64_t kappa;  // Used for the base-2^kappa representation (K)
	uint64_t n_limbs;  // (k+1)l or (k+1)l - 1 if l_a and l_b are different l,
	                   // at each limb we have a polynomial of degree N (Zn[X])
	double sigma;      // The standard deviation of the error distribution.
} GLWECtParams;

/**
 * @brief Creates a new set of bivGLWE parameters.
 *
 * @param N The polynomial maximum degree in X
 * @param k The number of Zn[X] polynomial in the secret key.
 * @param kappa The exponent for the base-2^kappa representation.
 * @param n_limbs (k+1)*l.
 * @param sigma The standard deviation of the error distribution.
 * @return GLWECtParams*
 */
GLWECtParams* new_glwe_ct_params(uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs, double sigma);

/**
 * @brief Deletes the set of bivGLWE parameters.
 *
 * @param params The bivGLWE parameters.
 */
void delete_glwe_ct_params(GLWECtParams* params);

#endif  // bivGLWE_CT_PARAMS_H
