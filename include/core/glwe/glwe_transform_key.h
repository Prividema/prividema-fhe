#ifndef bivGLWE_TRANSFORM_KEY_H
#define bivGLWE_TRANSFORM_KEY_H

#include "glwe_key.h"

/**
 * @brief Computes the secret key in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result bivGLWE secret key in the DFT domain.
 * @param sk The secret key out of the DFT domain.
 */
void transform_glwe_secret_key_not_dft_to_dft(const MODULE* module, GLWESecretKeyDFT* result_dft,
                                              const GLWESecretKey* sk);

#endif  //GLWE_TRANSFORM_KEY_H
