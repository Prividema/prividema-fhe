#ifndef GGSW_UTILS_H
#define GGSW_UTILS_H

#include "ggsw_ciphertext.h"
#include "glwegadget_ciphertext.h"

void check_ggsw(const MODULE* module, const GGSWCiphertext* ggsw, const GLWESecretKeyPrepared* sk_prep,
                const PolyUniv* expected, double max_err_length, double critical_err_length);

void check_glwegadget(const MODULE* module, const GLWEGadgetCiphertext* glwegad, const GLWESecretKeyPrepared* sk_prep,
                      const PolyUniv* expected, double max_err_length, double critical_err_length);
#endif
