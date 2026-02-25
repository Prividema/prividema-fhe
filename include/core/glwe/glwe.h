#ifndef GLWE_H
#define GLWE_H

#include "glwe_key.h"
#include "vec_znx_arithmetic_private.h"

// tnx element R[X] mod X^N+1 mod 1
typedef struct tnx_element
{
	uint64_t N;
	uint64_t kappa;
	uint64_t nb_limbs;
	double* coeffs;
} TNXElement;

//! GLWE PART (begin)
/**
 * @brief Computes Sum_j{0,k-1}[sk_j * a_j]
 *
 * @param module The module.
 * @param res The bivariate result.
 * @param ct The GLWE ciphertext
 * @param sk_dft The Secret key in DFT space.
 * @return int
 */
int add_mult(MODULE* module, GLWECtParams* params, PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft);

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
int glwe_secret_masking(MODULE* module, GLWECiphertext* ct, GLWESecretKeyDFT* sk_dft, PolyBiv* phase);

/**
 * @brief Computes -Sum_j{0,k-1}[sk_j * a_j]
 *
 * @param module The module.
 * @param res The bivariate result.
 * @param ct The GLWE ciphertext
 * @param sk_dft The Secret key in DFT space.
 * @return int
 */
int sub_mult(MODULE* module, GLWECtParams* params, PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft);

/**
 * @brief Demasks the bivariate phase (message + noise) and puts it in res.
 *
 * @param module
 * @param res The bivariate phase.
 * @param sk_dft The secret key in DFT space.
 * @param ct The bivGLWE ciphertext.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_demasking(MODULE* module, PolyBiv* res, GLWESecretKeyDFT* sk_dft, GLWECiphertext* ct);

// public key encrypt
void glwe_public_masking(const Core* core,    // all params of the library: is fft or
                                              // ntt, all N that are used
                         GLWECiphertext* ct,  // ciphertext
                         GLWEPublicKey* pk,   // public key
                         TNXElement* phase    // message + noise
);

// public key decrypt
void glwe_public_demasking(const Core* core,          // all params of the library: is fft or
                                                      // ntt, all N that are used
                           GLWECiphertext* ct,        // ciphertext
                           GLWESecretKeyDFT* sk_dft,  // secret key: vec of size k
                           TNXElement* phase          // message + noise
);

// addition 2 glwe
void glwe_addition(const Core* core, GLWECiphertext* ct_out, GLWECiphertext* ct_in1, GLWECiphertext* ct_in2);

//! GLWE IN DFT PART (begin)

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
int glwe_secret_masking_dft(MODULE* module, GLWECiphertextDFT* ct_dft, GLWESecretKeyDFT* sk_dft, PolyBivDFT* phase_dft);

/**
 * @brief Demasks the bivariate phase (message + noise) and puts it in res.
 *
 * @param res The bivariate phase.
 * @param sk_dft The secret key in DFT space.
 * @param ct The bivGLWE ciphertext.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_demasking_dft(MODULE* module, PolyBiv* res, GLWESecretKeyDFT* sk_dft, GLWECiphertextDFT* ct);

// public key encrypt
void glwe_public_masking_dft(const Core* core,    // all params of the library: is fft or
                                                  // ntt, all N that are used
                             GLWECiphertext* ct,  // ciphertext
                             GLWEPublicKey* pk,   // public key
                             TNXElement* phase    // message + noise
);

// public key decrypt
void glwe_public_demasking_dft(const Core* core,          // all params of the library: is fft or
                                                          // ntt, all N that are used
                               GLWECiphertext* ct,        // ciphertext
                               GLWESecretKeyDFT* sk_dft,  // secret key: vec of size k
                               TNXElement* phase          // message + noise
);

// addition 2 glwe
void glwe_addition_dft(const Core* core, GLWECiphertext* ct_out, GLWECiphertext* ct_in1, GLWECiphertext* ct_in2);

#endif  // GLWE_H