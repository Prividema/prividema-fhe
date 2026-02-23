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
 * @param params The GLWE parameters.
 * @param res The result bivariate phase.
 * @param phase The input phase.
 */
void add_error(MODULE* module, 
               GLWECtParams* params,
               PolyBiv* res,
               PolyBiv* phase);

/**
 * @brief Demasks the phase (message + noise).
 * 
 * @param params The GLWE parameters.
 * @param phase The result phase in Zn[X,Y]. 
 * @param sk_dft The secret key in DFT space.
 * @param ct The ciphertext.
 */
int glwe_secret_demasking_ggsw_lib(MODULE* module, 
                                   GLWECtParams* params,
                                   PolyBiv* phase,
                                   GGSWSecretKeyDFT* sk_dft,
                                   VecBiv* ct);

/**
 * @brief Masks the phase (message + noise) and puts it in res.
 * 
 * @param module The module stocking the degree N.
 * @param params The GLWE parameters.
 * @param res_ct The result bivariate ciphertext. 
 * @param sk_dft The secret key in DFT space.
 * @param phase message + noise.
 * 
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 */
int glwe_secret_masking_ggsw_lib(MODULE* module,
                                 GLWECtParams* params,
                                 VecBiv* res_ct,
                                 GGSWSecretKeyDFT* sk,
                                 PolyBiv* phase);

                
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
int compute_phase_ij(MODULE* module, GGSWCtParams* params, 
                     GGSWSecretKeyDFT* sk_dft, PolyUniv* msg_univ, PolyUnivDFT* msg_univ_dft,
                     PolyBiv* phase_biv, PolyUnivRnX* phase_univ_RnX, 
                     PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft, 
                     int64_t i, int64_t j);

/**
 * @brief Encrypts the message m into GGSW ciphertext res with parameters params.
 * 
 * @param params The encryption params
 * @param res The encrypted message
 * @param sk_dft The secret key
 * @param msg_univ The message
 * 
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int ggsw_secret_encrypt(MODULE* module, 
                        GGSWCtParams* params,
                        GGSWCiphertext* ct_ggsw,           
                        GGSWSecretKeyDFT* sk_dft,             
                        PolyUniv* msg_univ);

/**
 * @brief Computes the external product between a bivGLWE and a bivGGSW.
 * 
 * @param res The bivariate GLWE result ciphertext.
 * @param ct_glwe The bivariate GLWE input ciphertext.
 * @param ct_ggsw The bivariate GGSW input ciphertext. 
 */
int ggsw_external_product(MODULE* module, 
                           GLWECiphertext* res,  // result
                           GLWECiphertext* ct_glwe,  // GLWE ciphertext
                           GGSWCiphertext* ct_ggsw   // GGSW ciphertext
);
// TODO : the 8 functions below are not implemented

// GGSWPublicKey is a struct encapsulating everything regarding the public
// key.
/* Encrypts message m into GGSW ciphertext res with parameters params */
void ggsw_public_encrypt(GGSWCiphertext* res,           // result
                         GGSWPublicKey* pk,             // secret key
                         int64_t* m,                // message
                         GGSWCtParams* params  // parameters
);

// PartialGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_secret_encrypt(PartialGGSWCiphertext* res,       // result
                             GGSWSecretKey* sk,             // secret key
                             int64_t* m,                // message
                             GGSWCtParams* params  // parameters
);

/* Encrypts message m into halfGGSW ciphertext res with parameters params */
void halfggsw_public_encrypt(PartialGGSWCiphertext* res,       // result
                             GGSWPublicKey* pk,             // public key
                             int64_t* m,                // message
                             GGSWCtParams* params  // parameters
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void halfggsw_decrypt(int64_t* res,       // result
                      GGSWSecretKey* sk,      // secret key
                      PartialGGSWCiphertext* ct  // ciphertext
);




//! GGSW IN DFT PART (begin)   

/**
 * @brief Adds a bivariate error to the bivariate phase, returns in DFT space.
 * 
 * @param params The GLWE parameters.
 * @param res_dft The result bivariate phase in DFT space.
 * @param phase_dft The input phase.
 */
void add_error_dft(MODULE* module, 
                   GLWECtParams* params,
                   PolyBivDFT* res_dft,
                   PolyBivDFT* phase_dft);

/**
 * @brief Computes Sum_j{0,k-1}[sk_j * a_j].
 * 
 * @param module 
 * @param params 
 * @param res 
 * @param as 
 * @param sk_dft 
 * @return int 
 */
int add_mult_dft_ggsw(MODULE* module, GLWECtParams* params,
                      PolyBiv* res, VecBiv* as, GGSWSecretKeyDFT* sk_dft
);

/**
 * @brief Demasks the phase (message + noise) in DFT space and computes it out of DFT space.
 * 
 * @param params The GLWE parameters.
 * @param phase The result phase in Zn[X,Y]. 
 * @param sk_dft The secret key in DFT space.
 * @param ct The ciphertext.
 */
int glwe_secret_demasking_ggsw_lib_dft(MODULE* module,
                                       GLWECtParams* params,
                                       PolyBiv* phase,
                                       GGSWSecretKeyDFT* sk_dft,
                                       VecBivDFT* ct_dft);

/**
 * @brief Masks the phase (message + noise) in DFT space and puts it in res_ct.
 * 
 * @param params The GLWE parameters.
 * @param module The module stocking the degree N.
 * @param res_dft The result ciphertext in DFT space. 
 * @param sk_dft The secret key in DFT space.
 * @param phase_dft message + error.
 * 
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 */
int glwe_secret_masking_ggsw_lib_dft(MODULE* module,
                                     GLWECtParams* params,
                                     VecBivDFT* res_dft,
                                     GGSWSecretKeyDFT* sk_dft,
                                     PolyBivDFT* phase_dft);

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
int compute_phase_ij_dft(MODULE* module, GGSWCtParams* params_ggsw, 
                         GGSWSecretKeyDFT* sk_dft, PolyUniv* msg_univ, PolyUnivDFT* msg_univ_dft,
                         PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft, 
                         PolyBivDFT* phase_dft, PolyBiv* phase, PolyUnivRnX* phase_univ_RnX, 
                         int64_t i, int64_t j);

/**
 * @brief Encrypts the message m into GGSW ciphertext res with parameters params in DFT space.
 * 
 * @param params The encryption params
 * @param res_dft The encrypted message
 * @param sk_dft The secret key
 * @param msg_univ The message
 * 
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int ggsw_secret_encrypt_dft(MODULE* module,
                            GGSWCtParams* params,
                            GGSWCiphertextDFT* ct_ggsw_dft,
                            GGSWSecretKeyDFT* sk_dft,
                            PolyUniv* msg_univ);

/**
 * @brief Computes the external product between a bivGLWE and a biv GGSW.
 * 
 * @param res_dft The bivariate GLWE result ciphertext in DFT space.
 * @param ct_glwe_dft The bivariate GLWE input ciphertext in DFT space.
 * @param ct_ggsw_dft The bivariate GGSW input ciphertext in DFT space.
 */
int ggsw_external_product_dft(MODULE* module,
                               GLWECiphertextDFT* res_dft,  
                               GLWECiphertextDFT* ct_glwe_dft, 
                               GGSWCiphertextDFT* ct_ggsw_dft  
);

//! COMMON PART (begin)

#endif  // GGSW_H
