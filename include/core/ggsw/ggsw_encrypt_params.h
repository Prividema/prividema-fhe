#ifndef GGSW_ENCRYPT_PARAMS_H
#define GGSW_ENCRYPT_PARAMS_H

#include <stdint.h>

#include "glwe_encrypt_params.h"

typedef struct ggsw_encrypt_params {
  GLWEEncryptParams* params;
 //int64_t l;  // Decomposition size for the gadget
 //int64_t B;  // Decomposition basis for the gadget
  uint64_t k_tilde; // k_tilde=1 for RGSW
  uint64_t kappa_tilde;  // B= 2^-kappa
  uint64_t n_limbs_tilde; 
} GGSWEncryptParams;

typedef struct partialggsw_encrypt_params {
  GLWEEncryptParams* params;
  uint64_t kappa_tilde;  // B= 2^-kappa
  uint64_t l_tilde; 
} PartialGGSWEncryptParams;

#endif  // GGSW_ENCRYPT_PARAMS_H
