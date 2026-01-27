#ifndef GGSW_KEY_H
#define GGSW_KEY_H

#include "ggsw_ciphertext.h"

typedef struct ggsw_secret_key {
  uint64_t N;
  uint64_t k;
  int64_t* values;     // The key itself.
} GGSWSecretKey;

typedef struct ggsw_prep_secret_key {
  uint64_t N;
  uint64_t k;
  PolyUnivDFT** values;  // vec of size k, each element is prepared vec
  void* data;
} GGSWPreparedSK;

/**
 * @brief Creates a GGSW Secret key in DFT space
 * 
 * @param values The values of the secret key in DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of Zn[X] polynomial in the secret key.
 * @return GGSWPreparedSK* 
 */
GGSWPreparedSK* new_ggsw_secret_key_dft(PolyUnivDFT** values, uint64_t N, uint64_t k);

/**
 * @brief Delete the secret key that is in DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 */
void delete_ggsw_secret_key_dft(GGSWPreparedSK* sk_dft);

typedef struct ggsw_public_key {
  uint32_t size;       // number of ciphertexts in public key
  GGSWCiphertext* pk;  // Public key is multiple encryptions of 0
} GGSWPublicKey;

#endif  // GGSW_KEY_H
