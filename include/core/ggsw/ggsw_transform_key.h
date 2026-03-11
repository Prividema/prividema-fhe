#include "ggsw_key.h"

/**
 * @brief Computes the input secret key out of the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result The result GGSW secret key out of the DFT domain.
 * @param sk_dft A GGSW secret key in the DFT domain.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int transform_ggsw_secret_key_dft_to_not_dft(const MODULE* module, GGSWSecretKey* result, const GGSWSecretKeyDFT* sk_dft);

/**
 * @brief Transforms a GGSW secret key to a GLWE secret key.
 *
 * @param result The result GLWE secret key.
 * @param sk_ggsw A GGSW secret key
 */
void transform_ggsw_secret_key_to_glwe_secret_key(GLWESecretKey* result, const GGSWSecretKey* sk_ggsw);


/**
 * @brief Computes the secret key in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param result_dft The result GGSW secret key in the DFT domain.
 * @param sk The secret key out of the DFT domain.
 */
void transform_ggsw_secret_key_not_dft_to_dft(const MODULE* module, GGSWSecretKeyDFT* result_dft, const GGSWSecretKey* sk);

/**
 * @brief Transforms a GGSW secret key in the DFT domain to a GLWE secret key in the DFT domain.
 *
 * @param result_dft The result GLWE secret key in the DFT domain.
 * @param sk_ggsw_dft The GGSW secret key in the DFT domain.
 * @return GLWESecretKey*
 */
void transform_ggsw_secret_key_dft_to_glwe_secret_key_dft(GLWESecretKeyDFT* result_dft, const GGSWSecretKeyDFT* sk_ggsw_dft);
