#ifndef bivGGSW_H
#define bivGGSW_H

#include "ggsw_key.h"
#include "spqlios_alias.h"

#define WITH_Y0 // TODO : Do we have to keep this ?

// =============================================
// |                                           |
// |      			 iDFT Form       		   |
// |                                           |
// =============================================

/**
 * TODO: move to GLWE code
 * @brief Adds a bivariate error to a bivariate polynomial.
 *
 * @param module      Additionnal information for backend.
 * @param params_glwe The GLWE parameters.
 * @param res         The result bivariate phase.
 * @param pol         The input polynomial.
 * 
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int add_bivariate_error(const MODULE* module, const GLWECtParams* params, PolyBiv* result, const PolyBiv* pol);

/**
 * @brief Demasks the phase (message + noise).
 *
 * @param module      Additionnal information for backend.
 * @param params_glwe The GLWE parameters.
 * @param result      The result phase in \f$\mathbb{Z}_n[X,Y]\f$.
 * @param sk_dft      The secret key in the DFT domain.
 * @param glwe_vec    The ciphertext.
 * 
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int glwe_secret_demasking_ggsw_lib(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result, const GGSWSecretKeyDFT* sk_dft, const VecBiv* glwe_vec);

/**
 * @brief Masks the phase (message + noise) and puts it in result.
 *
 * @param module      Additionnal information for backend.
 * @param params_glwe The GLWE parameters.
 * @param result      The result bivariate ciphertext.
 * @param sk_dft      The secret key in the DFT domain.
 * @param phase       message + noise.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int glwe_secret_masking_ggsw_lib(const MODULE* module, const GLWECtParams* params_glwe, VecBiv* result, const GGSWSecretKeyDFT* sk_dft,
                                 const PolyBiv* phase);

/**
 * @brief Encrypts the message m into bivGGSW ciphertext res with parameters params.
 *
 * @param module      Additionnal information for backend.
 * @param params_ggsw The encryption params.
 * @param result      The encrypted message.
 * @param sk_dft      The secret key.
 * @param m_univ      The message.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int ggsw_secret_encrypt(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertext* result, const GGSWSecretKeyDFT* sk_dft,
                        const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a bivGGSW.
 *
 * @param module Additionnal information for backend.
 * @param result The bivariate GLWE result ciphertext.
 * @param glwe   The bivariate GLWE input ciphertext.
 * @param ggsw   The bivariate GGSW input ciphertext.
 * 
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int ggsw_external_product(const MODULE* module, GLWECiphertext* result,
                          const GLWECiphertext* glwe, const GGSWCiphertext* ggsw);

// =============================================
// |                                           |
// |      			 DFT Form       		   |
// |                                           |
// =============================================

/**
 * @brief Demasks the phase (message + noise) in the DFT domain and computes it out of the DFT domain.
 *
 * @param module       Additionnal information for backend.
 * @param params_glwe  The GLWE parameters.
 * @param result       The result phase in \f$\mathbb{Z}_n[X,Y]\f$.
 * @param sk_dft       The secret key in the DFT domain.
 * @param glwe_vec_dft The ciphertext.
 * 
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int glwe_secret_demasking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result, const GGSWSecretKeyDFT* sk_dft,
                                       const VecBivDFT* glwe_vec_dft);

/**
 * @brief Masks the phase (message + noise) in the DFT domain and puts it in result.
 *
 * @param module      Additionnal information for backend.
 * @param params_glwe The GLWE parameters.
 * @param result_dft  The result ciphertext in the DFT domain.
 * @param sk_dft      The secret key in the DFT domain.
 * @param phase_dft   message + error.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int glwe_secret_masking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params_glwe, VecBivDFT* result_dft, const GGSWSecretKeyDFT* sk_dft,
                                     const PolyBivDFT* phase_dft);


/**
 * @brief Encrypts the message m into bivGGSW ciphertext res with parameters params in the DFT domain.
 *
 * @param module      Additionnal information for backend.
 * @param params_ggsw The encryption params.
 * @param res_dft     The encrypted message.
 * @param sk_dft      The secret key.
 * @param m_univ      The message.
 *
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int ggsw_secret_encrypt_dft(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertextDFT* result_dft,
                            const GGSWSecretKeyDFT* sk_dft, const PolyUniv* m_univ);

/**
 * @brief Computes the external product between a bivGLWE and a biv bivGGSW.
 *
 * @param module   Additionnal information for backend.
 * @param res_dft  The bivariate GLWE result ciphertext in the DFT domain.
 * @param glwe_dft The bivariate GLWE input ciphertext in the DFT domain.
 * @param ggsw_dft The bivariate GGSW input ciphertext in the DFT domain.
 * 
 * @retval • `-1` if an error occurs.
 * @retval • `0` otherwise otherwise.
 */
int ggsw_external_product_dft(const MODULE* module,
                              GLWECiphertextDFT* result_dft,     
                              const GLWECiphertextDFT* glwe_dft, 
                              const GGSWCiphertextDFT* ggsw_dft);

#endif  // bivGGSW_H
