#ifndef GGSW_CT_PARAMS_H
#define GGSW_CT_PARAMS_H

#include <stdint.h>

#include "glwe_ct_params.h"

typedef struct ggsw_ct_params {
  GLWECtParams* params;
  uint64_t k_tilde; // k_tilde=1 for RGSW (by default k=k_tilde=1)
  uint64_t kappa_tilde;  // B= 2^-kappa
  uint64_t n_limbs_tilde; 
} GGSWCtParams;

typedef struct partialggsw_ct_params {
  GLWECtParams* params;
  uint64_t kappa_tilde;  // B= 2^-kappa
  uint64_t l_tilde; 
} PartialGGSWCtParams;

#endif  // GGSW_CT_PARAMS_H
