#ifndef GLWE_ENCRYPT_PARAMS_H
#define GLWE_ENCRYPT_PARAMS_H

#include <stdint.h>

typedef struct glwe_encrypt_params {
  int64_t N;     // Polynomial degree
  int64_t k;     // Number of a terms
  int64_t base;  // Base used for the base-2^k representation
  int64_t prec;  // Decomposition size for base-2^k representation
} GLWEEncryptParams;

#endif  // GLWE_ENCRYPT_PARAMS_H
