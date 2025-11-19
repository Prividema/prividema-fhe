#ifndef GGSW_KEY_H
#define GGSW_KEY_H

#include "ggsw_ciphertext.h"


typedef struct ggsw_secret_key {
  uint32_t size;        // Key size, usually k.N
  int64_t* values;  // The key itself.
} GGSWSecretKey;

typedef struct ggsw_public_key {
  uint32_t size;       // number of ciphertexts in public key
  GGSWCiphertext* pk;  // Public key is multiple encryptions of 0
} GGSWPublicKey;

#endif  // GGSW_KEY_H
