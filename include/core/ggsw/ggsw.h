#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_ciphertext.h"
#include "ggsw_key.h"
#include "spqlios_alias.h"

#define WITH_Y0

//! bivGGSW PART (begin)

/**
 * @brief Adds a bivariate error to the bivariate phase.
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param res The result bivariate phase.
 * @param phase The input bivariate phase.
 */
int add_error(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result, const PolyBiv* phase);

/**
 * @brief Demasks the phase (message + noise).
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param result The result bivariate phase.
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe_vec The bivGLWE ciphertext.
 */
int glwe_secret_demasking_ggsw_lib(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result, const GGSWSecretKeyDFT* sk_dft, const VecBiv* glwe_vec);

/**
 * @brief Masks the phase (message + noise) and puts it in result.
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase message + noise.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise otherwise.
 */
int glwe_secret_masking_ggsw_lib(const MODULE* module, const GLWECtParams* params_glwe, VecBiv* result, const GGSWSecretKeyDFT* sk_dft,
                                 const PolyBiv* phase);

/**
 * @brief Computes the base-2kappa decomposition of the univariate phase : -m * sk_j / 2^{kappa_tilde * i} + err, if j < k
                                                                            m / 2^{kappa_tilde * i} + err       , if j = k
 *
 * @param module Additionnal information for backend.
 * @param params_ggsw The bivGGSW parameters.
 * @param result The bivariate result phase.
 * @param phase_univ_RnX The phase in univariate space.
 * @param m_skj_univ The product m * sk_j, for j < k.
 * @param m_skj_univ_dft The product m * sk_j in the DFT domain, for j < k.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ_dft The univariate message in the DFT domain.
 * @param m_univ The univariate message.
 * @param i The i-ème block of HalfGGSW = [GLWE(-m*sk_0 / 2^{kappa_tilde * i}), ... , bivGLWE(-m*sk_(k-1) / 2^{kappa_tilde * i}), bivGLWE(m / 2^{kappa_tilde * i})]
 * @param j The index of the element of the secret key.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise otherwise.
 */
int compute_phase_ij(const MODULE* module, const GGSWCtParams* params_ggsw, PolyBiv* result, 
					 PolyUnivRnX* phase_univ_RnX, PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft, 
					 const GGSWSecretKeyDFT* sk_dft, const PolyUnivDFT* m_univ_dft, const PolyUniv* m_univ, 
                     int64_t i, int64_t j);

/**
 * @brief Encrypts the message m into bivGGSW ciphertext res with parameters params.
 *
 * @param module Additionnal information for backend.
 * @param params_ggsw The encryption params.
 * @param result The bivGGSW(m) encrypted message.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ The univariate message.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise otherwise.
 */
int ggsw_secret_encrypt(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertext* result, const GGSWSecretKeyDFT* sk_dft,
                        const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a bivGGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The bivGLWE result ciphertext.
 * @param glwe The bivGLWE input ciphertext.
 * @param ggsw The bivGGSW input ciphertext.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise otherwise.
 */
int ggsw_external_product(const MODULE* module,
                          GLWECiphertext* result,      // result
                          const GLWECiphertext* glwe,  // bivGLWE ciphertext
                          const GGSWCiphertext* ggsw   // bivGGSW ciphertext
);

// PartialGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_secret_encrypt(PartialGGSWCiphertext* result,  // result
                             const GGSWSecretKey* sk,           // secret key
                             const int64_t* m,                  // message
                             const GGSWCtParams* params         // parameters
);


//! bivGGSW IN DFT PART (begin)

/**
 * @brief Demasks the phase (message + noise) in the DFT domain and computes it out of the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param result result phase in Zn[X,Y].
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe_vec_dft The ciphertext.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int glwe_secret_demasking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result, const GGSWSecretKeyDFT* sk_dft,
                                       const VecBivDFT* glwe_vec_dft);

/**
 * @brief Masks the phase (message + noise) in the DFT domain and puts it in result.
 *
 * @param module Additionnal information for backend.
 * @param params_glwe The bivGLWE parameters.
 * @param result_dft The result ciphertext in the DFT domain.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase_dft message + error.
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise otherwise.
 */
int glwe_secret_masking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params_glwe, VecBivDFT* result_dft, const GGSWSecretKeyDFT* sk_dft,
                                     const PolyBivDFT* phase_dft);

/**
 * @brief Computes the base-2kappa decomposition of the phase : -m * sk_j / 2^{kappa_tilde^i} + err, if j < k
 *                                                                      m / 2^{kappa_tilde^i} + err, if j = k
 *
 * @param module Additionnal information for backend.
 * @param params_ggsw The bivGGSW parameters.
 * @param sk_dft The secret key in the DFT domain.
 * @param m_univ The input message.
 * @param m_univ_dft The input message in the DFT domain.
 * @param phase_dft The phase in the DFT domain.
 * @param phase The phase.
 * @param phase_univ_RnX
 * @param m_skj_univ
 * @param m_skj_univ_dft
 * @param i
 * @param j
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int compute_phase_ij_dft(const MODULE* module, const GGSWCtParams* params_ggsw, PolyBivDFT* result_dft, 
						 PolyBiv* phase, PolyUnivRnX* phase_univ_RnX, PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft, 
						 const GGSWSecretKeyDFT* sk_dft, const PolyUnivDFT* m_univ_dft, const PolyUniv* m_univ,
						 int64_t i, int64_t j);

/**
 * @brief Encrypts the message m into bivGGSW ciphertext res with parameters params in the DFT domain.
 *
 * @param module Additionnal information for backend.
 * @param params_ggsw The encryption params
 * @param res_dft The encrypted message
 * @param sk_dft The secret key
 * @param m_univ The message
 *
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise otherwise.
 */
int ggsw_secret_encrypt_dft(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertextDFT* result_dft,
                            const GGSWSecretKeyDFT* sk_dft, const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a biv bivGGSW.
 *
 * @param module Additionnal information for backend.
 * @param res_dft The bivariate bivGLWE result ciphertext in the DFT domain.
 * @param glwe_dft The bivariate bivGLWE input ciphertext in the DFT domain.
 * @param ggsw_dft The bivariate bivGGSW input ciphertext in the DFT domain.
 */
int ggsw_external_product_dft(const MODULE* module,
                              GLWECiphertextDFT* result_dft,     
                              const GLWECiphertextDFT* glwe_dft, 
                              const GGSWCiphertextDFT* ggsw_dft);
//! COMMON PART (begin)

#endif  // bivGGSW_H
