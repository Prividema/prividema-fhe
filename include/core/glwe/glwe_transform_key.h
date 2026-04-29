#ifndef bivGLWE_TRANSFORM_KEY_H
#define bivGLWE_TRANSFORM_KEY_H

#include "glwe_key.h"

/**
 * @brief Transforms a secret key from coefficient space to DFT
 *
 * @param module Additionnal information for backend.
 * @param result_dft The resulting GLWE secret key in the DFT domain.
 * @param sk The secret key in coefficient space
 */
void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft,
                                              const GLWESecretKey* sk);

#endif  //GLWE_TRANSFORM_KEY_H
