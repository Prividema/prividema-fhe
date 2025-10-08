#ifndef HALFGGSW_CIPHERTEXT_H
#define HALFGGSW_CIPHERTEXT_H

#include "ggsw_encrypt_params.h"
#include "glwe_ciphertext.h"

typedef struct halfggsw_ciphertext {
  GGSWEncryptParams* params;
  GLWECiphertext* values;  // vector of RLWE
} HalfGGSWCiphertext;

#endif  // HALFGGSW_CIPHERTEXT_H
