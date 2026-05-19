#ifndef GGSW_UTILS_H
#define GGSW_UTILS_H

#include "ggsw_ciphertext.h"

void check_ggsw(const MODULE* module, const GLWEParams* params_glwe, const GGSWParams* params_ggsw,
                const GGSWCiphertext* ggsw, const GLWESecretKeyPrepared* sk_prep, const PolyUniv* expected,
                double max_err_length, double critical_err_length);

#endif
