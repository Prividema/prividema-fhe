#ifndef GGSW_UTILS_H
#define GGSW_UTILS_H

#include "ggsw_ciphertext.h"

/**
 * @brief Test utility function to check that a GGSW is a valid and "close enough" ciphertext of a certain plaintext
 *
 * @param module        The backend module
 * @param ggsw          The GGSW ciphertext under test
 * @param sk_prep       The secret key used for encryption, in prepared form
 * @param expected      The expected plaintext
 * @param max_err       Error tolerance as in pvda_assert_polynomial_distance
 * @param critical_err  Error tolerance as in pvda_assert_polynomial_distance
 *
 */
void check_ggsw(const PvdaBackend* module, const GGSWCiphertext* ggsw, const GLWESecretKeyPrepared* sk_prep,
                const PolyUniv* expected, double max_err, double critical_err);

#endif
