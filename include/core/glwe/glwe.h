#ifndef GLWE_H
#define GLWE_H

#include "glwe_key.h"

// tnx element R[X] mod X^N+1 mod 1
typedef struct tnx_element {
  uint64_t N;
  uint64_t kappa;
  uint64_t nb_limbs;
  double* coeffs;
} TNXElement;

// secret key encrypt
int glwe_secret_masking(GLWECiphertext* ct,  // ciphertext
                      GLWEPreparedSK* sk_dft,   // secret key: vec of size k
                      PolyBiv* phase    // message + noise
);

// secret key decrypt (compute the phase)
int glwe_secret_demasking(GLWECtParams* enc_params,
                          TNXElement* phase,  
                          GLWEPreparedSK* sk_dft, 
                          GLWECiphertext* ct 
);

// add noise message

// public key encrypt
void glwe_public_masking(const Core* core,  // all params of the library: is fft or
                                         // ntt, all N that are used
                      GLWECiphertext* ct,  // ciphertext
                      GLWEPublicKey* pk,   // public key
                      TNXElement* phase    // message + noise
);

// public key decrypt
void glwe_public_demasking(const Core* core,    // all params of the library: is fft or
                                         // ntt, all N that are used
                    GLWECiphertext* ct,  // ciphertext
                    GLWEPreparedSK* sk_dft,   // secret key: vec of size k
                    TNXElement* phase    // message + noise
);

// addition 2 glwe
void glwe_addition(const Core* core, GLWECiphertext* ct_out,
                   GLWECiphertext* ct_in1, GLWECiphertext* ct_in2);

#endif  // GLWE_H