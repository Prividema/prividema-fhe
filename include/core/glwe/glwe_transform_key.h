#ifndef GLWE_TRANSFORM_KEY_H
#define GLWE_TRANSFORM_KEY_H

#include "glwe_key.h"

//! GLWE KEY PART (begin)

/**
 * @brief Computes the secret key out of the DFT domain.
 *
 * @param sk_dft The secret key in the DFT domain.
 * @return GLWESecretKey*
 */
GLWESecretKey* transform_glwe_secret_key_dft_to_not_dft(const MODULE* module, const GLWESecretKeyDFT* sk_dft);

/**
 * @brief Computes the values of the secret key out of the DFT domain.
 *
 * @param values_dft The values of the secret key in the DFT domain.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUniv**
 */
PolyUniv** transform_glwe_secret_key_values_dft_to_not_dft(const MODULE* module, const PolyUnivDFT** values_dft, uint64_t k);


//! GLWE IN DFT SPACE PART (begin)

/**
 * @brief Computes the secret key in the DFT domain.
 *
 * @param sk The secret key out of the DFT domain.
 * @return GLWESecretKey*
 */
GLWESecretKeyDFT* transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, const GLWESecretKey* sk);

/**
 * @brief Computes the values of the secret key in the DFT domain.
 *
 * @param values The values of the secret key out of the DFT domain.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return PolyUnivDFT**
 */
PolyUnivDFT** transform_glwe_secret_key_values_not_dft_to_dft(const MODULE* module, const PolyUniv** values, uint64_t k);

#endif //GLWE_TRANSFORM_KEY_H