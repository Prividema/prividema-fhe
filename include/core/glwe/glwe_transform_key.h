#ifndef GLWE_TRANSFORM_KEY_H
#define GLWE_TRANSFORM_KEY_H

#include "glwe_key.h"

//! GLWE KEY PART (begin)

/**
 * @brief Computes the input secret key out of the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result The result GLWE secret key out of the DFT domain.
 * @param sk_dft A GLWE secret key in the DFT domain.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int transform_glwe_secret_key_dft_to_not_dft(const MODULE* module, GLWESecretKey* result, const GLWESecretKeyDFT* sk_dft);


//! GLWE IN DFT SPACE PART (begin)

/**
 * @brief Computes the secret key in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result GLWE secret key in the DFT domain.
 * @param sk The secret key out of the DFT domain.
 */
void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft, const GLWESecretKey* sk);

#endif //GLWE_TRANSFORM_KEY_H