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

/**
 * @brief Masks the phase (message + noise) and puts it in res.
 *  
 * @param ct The result bivariate ciphertext. 
 * @param sk_dft The secret key in DFT space.
 * @param phase message + noise.
 * 
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_masking(GLWECiphertext* ct,  // ciphertext
                        GLWEPreparedSK* sk_dft,   // secret key: vec of size k
                        PolyBiv* phase    // message + noise
);

/**
 * @brief Decrypts the phase (message + noise) and puts it in phase.
 * 
 * @param enc_params The GLWE parameters.
 * @param phase The phase in Zn[X,Y]. 
 * @param sk_dft The secret key in DFT space.
 * @param ct The ciphertext.
 * 
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_demasking(GLWECtParams* enc_params,
                          PolyBiv* phase,  
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