#ifndef GLWE_CT_PARAMS_H
#define GLWE_CT_PARAMS_H

#include <stdint.h>
#include <stdlib.h>

typedef struct glwe_ct_params {
  uint64_t N;         // Polynomial degree
  uint64_t k;         // Number of a terms, k=1 for RLWE
  uint64_t kappa;     // Used for the base-2^kappa representation (K)
  // ct(a,b)
  uint64_t n_limbs;   // (k+1)l or (k+1)l - 1 if l_a and l_b are different l, 
                      // at each limb we have a polynomial of degree N (ZnX)
} GLWECtParams;

int new_glwe_ct_params(GLWECtParams* res, uint64_t N, uint64_t k, uint64_t kappa, uint64_t n_limbs);
void delete_glwe_ct_params(GLWECtParams* params);
int64_t glwe_size(GLWECtParams* params);
int64_t poly_biv_size(GLWECtParams* params);
#endif  // GLWE_CT_PARAMS_H
