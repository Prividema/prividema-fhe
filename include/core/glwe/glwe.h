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
 * @param module Additionnal information for backend.
 * @param result The bivariate polynomial result.
 * @param glwe The GLWE ciphertext
 * @param sk_dft The Secret key in the DFT domain.
 * @return int 
 */
int add_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* result, VecBiv* glwe, GLWESecretKeyDFT* sk_dft);

/**
 * @brief Masks the phase (message + noise) and puts it in result.
 *
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase message + noise.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_masking(const MODULE* module, GLWECiphertext* result, GLWESecretKeyDFT* sk_dft, PolyBiv* phase);

/**
 * @brief Computes -Sum_j{0,k-1}[sk_j * a_j]
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate result.
 * @param glwe The GLWE ciphertext
 * @param sk_dft The Secret key in the DFT domain.
 * @return int
 */
int sub_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* result, VecBiv* glwe_vec, GLWESecretKeyDFT* sk_dft);

/**
 * @brief Demasks the bivariate phase (message + noise) and puts it in result.
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate phase.
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe The bivGLWE ciphertext.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_demasking(const MODULE* module, PolyBiv* result, GLWESecretKeyDFT* sk_dft, GLWECiphertext* glwe);

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


//! GLWE IN DFT PART (begin)

/**
 * @brief Masks the phase (message + noise) in the DFT domain.
 *
 * @param result_dft The result bivariate GLWE ciphertext in the DFT domain..
 * @param sk_dft The secret key in the DFT domain.
 * @param phase_dft message + noise in the DFT domain.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_masking_dft(const MODULE* module, GLWECiphertextDFT* result_dft, GLWESecretKeyDFT* sk_dft, PolyBivDFT* phase_dft);

/**
 * @brief Demasks the bivariate phase (message + noise) and puts it in result.
 *
 * @param result The bivariate phase.
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe_dft The bivGLWE ciphertext in the DFT domain.
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int glwe_secret_demasking_dft(const MODULE* module, PolyBiv* result, GLWESecretKeyDFT* sk_dft, GLWECiphertextDFT* ct);

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

#endif  // GLWE_H