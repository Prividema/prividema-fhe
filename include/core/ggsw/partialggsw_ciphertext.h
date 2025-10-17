#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include "ggsw_encrypt_params.h"
#include "glwe_ciphertext.h"

typedef struct partialggsw_ciphertext {
  PartialGGSWEncryptParams* params;
  int64_t* ct;  // matrix of size n_limbs x l_tilde
} PartialGGSWCiphertext;

#endif  // PARTIALGGSW_CIPHERTEXT_H
