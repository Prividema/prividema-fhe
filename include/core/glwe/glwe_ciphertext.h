#ifndef GLWE_CIPHERTEXT_H
#define GLWE_CIPHERTEXT_H

#include <stdint.h>

// TODO : To be Defined (spqlios).
typedef struct b {
} Base2kIntegerPoly;

typedef struct glwe_ciphertext {
  uint64_t N;
  uint64_t k; // k=1 for RLWE
  uint64_t kappa;  // base_2k
  uint64_t n_limbs; //(k+1)l or (k+1)l - 1 if l_a and l_b are differe,t l, at each n_limbs we have a polynomial of degree N
  // l_a = floor((n_limbs+1)/(k+1))
  // l_b= n_limbs- kl_a
  // log2alpha = l_b*kappa 
  int64_t* ct;  //  corresponds to (a_0, ..., a_{k-1}, b=a_k)
  // ct[p]= a_i[j], i= p%(k+1) in [0, l_a[ and j= floor(p/k+1) in [O,k[
} GLWECiphertext;

#endif  // GLWE_CIPHERTEXT_H
