#ifndef GLWE_ENCRYPT_PARAMS_H
#define GLWE_ENCRYPT_PARAMS_H

#include <stdint.h>

typedef struct glwe_encrypt_params {
  uint64_t N;     // Polynomial degree
  uint64_t k;     // Number of a terms, k=1 for RLWE
  uint64_t kappa;  // Used for the base-2^kappa representation
  uint64_t n_limbs; //(k+1)l or (k+1)l - 1 if l_a and l_b are different l, at each n_limbs we have a polynomial of degree N
} GLWEEncryptParams;

#endif  // GLWE_ENCRYPT_PARAMS_H
