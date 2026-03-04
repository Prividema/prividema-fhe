#ifndef GGSW_H
#define GGSW_H

#include "ggsw_ciphertext.h"
#include "ggsw_key.h"
#include "spqlios_alias.h"

#define WITH_Y0

//! GGSW PART (begin)

/**
 * @brief Adds a bivariate error to the bivariate phase.
 *
 * @param module TODO
 * @param params The GLWE parameters.
 * @param res The result bivariate phase.
 * @param phase The input phase.
 */
int add_error(const MODULE* module, const GLWECtParams* params, PolyBiv* result, const PolyBiv* phase);

/**
 * @brief Demasks the phase (message + noise).
 *
 * @param module TODO
 * @param params The GLWE parameters.
 * @param result The result phase in Zn[X,Y].
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe_vec The ciphertext.
 */
int glwe_secret_demasking_ggsw_lib(const MODULE* module, const GLWECtParams* params, PolyBiv* result, const GGSWSecretKeyDFT* sk_dft, const VecBiv* glwe_vec);

/**
 * @brief Masks the phase (message + noise) and puts it in res_ct.
 *
 * @param module The module stocking the degree N.
 * @param params The GLWE parameters.
 * @param result The result bivariate ciphertext.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase message + noise.
 *
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 */
int glwe_secret_masking_ggsw_lib(const MODULE* module, const GLWECtParams* params, VecBiv* result, const GGSWSecretKeyDFT* sk_dft,
                                 const PolyBiv* phase);

/**
 * @brief Computes the base-2kappa decomposition of the phase : -m * sk_j / 2^{kappa_tilde*(i+1)} + err, if j < k
                                                                        m / 2^{kappa_tilde*(i+1)} + err, if j = k
 *
 * @param module The module stocking the degree N.
 * @param params_ggsw The GGSW parameters.
 * @param params_glwe
 * @param ct_ggsw
 * @param sk_dft
 * @param msg_univ_dft
 * @param phase_biv
 * @param phase_univ_RnX
 * @param m_skj_univ
 * @param m_skj_univ_dft
 * @param i
 * @param j
 * @return int
 */
int compute_phase_ij(const MODULE* module, const GGSWCtParams* params_ggsw, const GGSWSecretKeyDFT* sk_dft, const PolyUniv* msg_univ,
                     const PolyUnivDFT* msg_univ_dft, PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft, 
                     PolyBiv* result, PolyUnivRnX* phase_univ_RnX, int64_t i, int64_t j);

/**
 * @brief Encrypts the message m into GGSW ciphertext res with parameters params.
 *
 * @param module TODO
 * @param params_ggsw The encryption params
 * @param result The encrypted message
 * @param sk_dft The secret key
 * @param msg_univ The message
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int ggsw_secret_encrypt(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertext* result, const GGSWSecretKeyDFT* sk_dft,
                        const PolyUniv* msg_univ);

/**
 * @brief Computes the external product between a bivGLWE and a bivGGSW.
 *
 * @param module TODO
 * @param result The bivariate GLWE result ciphertext.
 * @param glwe The bivariate GLWE input ciphertext.
 * @param ggsw The bivariate GGSW input ciphertext.
 */
int ggsw_external_product(const MODULE* module,
                          GLWECiphertext* result,      // result
                          const GLWECiphertext* glwe,  // GLWE ciphertext
                          const GGSWCiphertext* ggsw   // GGSW ciphertext
);

// The three functions below are not implemented

// GGSWPublicKey is a struct encapsulating everything regarding the public
// key.
/* Encrypts message m into GGSW ciphertext res with parameters params */
void ggsw_public_encrypt(GGSWCiphertext* result,  // result
                         const GGSWPublicKey* pk,    // secret key
                         const int64_t* m,           // message
                         const GGSWCtParams* params  // parameters
);

// PartialGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_secret_encrypt(PartialGGSWCiphertext* result,  // result
                             const GGSWSecretKey* sk,           // secret key
                             const int64_t* m,                  // message
                             const GGSWCtParams* params         // parameters
);

/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_public_encrypt(PartialGGSWCiphertext* result,  // result
                             const GGSWPublicKey* pk,           // public key
                             const int64_t* m,                  // message
                             const GGSWCtParams* params         // parameters
);


//! GGSW IN DFT PART (begin)

/**
 * @brief Demasks the phase (message + noise) in the DFT domain and computes it out of DFT space.
 *
 * @param module TODO
 * @param params The GLWE parameters.
 * @param result result phase in Zn[X,Y].
 * @param sk_dft The secret key in the DFT domain.
 * @param glwe_vec_dft The ciphertext.
 */
int glwe_secret_demasking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params, PolyBiv* result, const GGSWSecretKeyDFT* sk_dft,
                                       const VecBivDFT* glwe_vec_dft);

/**
 * @brief Masks the phase (message + noise) in the DFT domain and puts it in res_ct.
 *
 * @param module The module stocking the degree N.
 * @param params The GLWE parameters.
 * @param result_dft The result ciphertext in the DFT domain.
 * @param sk_dft The secret key in the DFT domain.
 * @param phase_dft message + error.
 *
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 */
int glwe_secret_masking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params, VecBivDFT* result_dft, const GGSWSecretKeyDFT* sk_dft,
                                     const PolyBivDFT* phase_dft);

/**
 * @brief Computes the base-2kappa decomposition of the phase : -m * sk_j / 2^{kappa_tilde*(i+1)} + err, if j < k
 *                                                                      m / 2^{kappa_tilde*(i+1)} + err, if j = k
 *
 * @param module
 * @param params_ggsw
 * @param ct_ggsw
 * @param sk_dft
 * @param msg_univ
 * @param msg_univ_dft
 * @param phase_dft
 * @param phase
 * @param phase_univ_RnX
 * @param m_skj_univ
 * @param m_skj_univ_dft
 * @param i
 * @param j
 * @return int
 */
int compute_phase_ij_dft(const MODULE* module, const GGSWCtParams* params_ggsw, const GGSWSecretKeyDFT* sk_dft, const PolyUniv* msg_univ,
                         const PolyUnivDFT* msg_univ_dft, PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft,
                         PolyBivDFT* phase_dft, PolyBiv* result, PolyUnivRnX* phase_univ_RnX, int64_t i, int64_t j);

/**
 * @brief Encrypts the message m into GGSW ciphertext res with parameters params in the DFT domain.
 *
 * @param module TODO
 * @param params The encryption params
 * @param res_dft The encrypted message
 * @param sk_dft The secret key
 * @param msg_univ The message
 *
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int ggsw_secret_encrypt_dft(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertextDFT* result_dft,
                            const GGSWSecretKeyDFT* sk_dft, const PolyUniv* msg_univ);

/**
 * @brief Computes the external product between a bivGLWE and a biv GGSW.
 *
 * @param module TODO
 * @param res_dft The bivariate GLWE result ciphertext in the DFT domain.
 * @param ct_glwe_dft The bivariate GLWE input ciphertext in the DFT domain.
 * @param ct_ggsw_dft The bivariate GGSW input ciphertext in the DFT domain.
 */
int ggsw_external_product_dft(const MODULE* module,
                              GLWECiphertextDFT* result_dft,     
                              const GLWECiphertextDFT* glwe_dft, 
                              const GGSWCiphertextDFT* ggsw_dft);
//! COMMON PART (begin)

#endif  // GGSW_H
