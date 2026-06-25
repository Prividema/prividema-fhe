#ifndef bivGLWE_TRANSFORM_KEY_H
#define bivGLWE_TRANSFORM_KEY_H

#include "glwe_key.h"

/**
 *
 * @file glwe_transform_key.h
 *
 * Operations to transform a GLWE/GLWEGadget/GGSW secret key (only coef to DFT is needed for now)
 */

/**
 * @brief Transforms a secret key from coefficient space to DFT
 *
 * @param module Additionnal information for backend.
 * @param result_dft The resulting GLWE secret key in the DFT domain.
 * @param sk The secret key in coefficient space
 */
void glwe_sk_prepare(const PVDA_MODULE* module, GLWESecretKeyPrepared* result_dft, const GLWESecretKey* sk);

#endif  //GLWE_TRANSFORM_KEY_H
