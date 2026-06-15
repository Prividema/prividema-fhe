#ifndef GLWEGAD_UTILS
#define GLWEGAD_UTILS

#include "ggsw_params.h"
#include "glwe_params.h"
#include "glwegadget_ciphertext.h"
#include "spqlios_alias.h"

/**
 * @brief Test utility function to check that a GLWEGadget is a valid and "close enough" ciphertext of a certain plaintext
 *
 * @param module        The backend module
 * @param glwegad       The GGSW ciphertext under test
 * @param sk_prep       The secret key used for encryption, in prepared form
 * @param expected      The expected plaintext
 * @param max_err       Error tolerance as in pvda_assert_polynomial_distance
 * @param critical_err  Error tolerance as in pvda_assert_polynomial_distance
 *
 */
void check_glwegadget(const MODULE* module, const GLWEGadgetCiphertext* glwegad, const GLWESecretKeyPrepared* sk_prep,
                      const PolyUniv* expected, double max_err, double critical_err);

/**
 * @brief Very rought estimate the number of bits that are not noisy after a half product
 *
 * Only takes into consideration the number of precision bits in both inputs.
 * Thus, it will not give a strictly correct result when there is any kind of noise in the
 * ciphertext, which is always the case except in test dummy cases
 *
 * @param params_glwe       The parameters of the GLWE being multiplied
 * @param params_glwegadget The parameters of the GLWEGadget being multiplied
 *
 */
int64_t info_bits_half_prod(const GLWEParams* params_glwe, const GLWEGadgetParams* params_glwegadget);

#endif
