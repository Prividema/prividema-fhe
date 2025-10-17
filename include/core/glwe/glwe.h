#ifndef GLWE_H
#define GLWE_H

#include "glwe_key.h"

// tnx element R[X] mod X^N+1 mod 1
typedef struct tnx_element {
  uint64_t N;
  uint64_t kappa;
  uint64_t l;
  int64_t* coeffs;
  void* data;
} TNXElement;

// secret key encrypt
void glwe_encrypt_priv(const Core* core,    // all params of the library: is fft
                                            // or ntt, all N that are used
                       GLWECiphertext* ct,  // ciphertext
                       GLWEPreparedSK* s,   // secret key: vec of size k
                       TNXElement* phase    // message + noise
);

// secret key decrypt (compute the phase)
void glwe_phase_priv(const Core* core,   // all params of the library: is fft or
                                         // ntt, all N that are used
                     TNXElement* phase,  // decrypted phase
                     GLWEPreparedSK* s,  // secret key
                     GLWECiphertext* ct  // ciphertext
);

// public key encrypt
void glwe_encrypt_pub(const Core* core,  // all params of the library: is fft or
                                         // ntt, all N that are used
                      GLWECiphertext* ct,  // ciphertext
                      GLWEPublicKey* pk,   // public key
                      TNXElement* phase    // message + noise
);

// public key decrypt
void glwe_phase_pub(const Core* core,    // all params of the library: is fft or
                                         // ntt, all N that are used
                    GLWECiphertext* ct,  // ciphertext
                    GLWEPreparedSK* s,   // secret key: vec of size k
                    TNXElement* phase    // message + noise
);

// addition 2 glwe
void glwe_addition(const Core* core, GLWECiphertext* ct_out,
                   GLWECiphertext* ct_in1, GLWECiphertext* ct_in2);

#endif  // GLWE_H
