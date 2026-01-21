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
 * @param enc_params The GLWE parameters.
 * @param res The result bivariate phase.
 * @param phase The input phase.
 */
void add_error(GLWECtParams* enc_params,
               PolyBiv* res,
               PolyBiv* phase);

/**
 * @brief Decrypts the phase (message + noise) and puts it in phase.
 * 
 * @param enc_params The GLWE parameters.
 * @param phase The phase in Rn[X]. 
 * @param sk_dft The secret key in DFT space.
 * @param ct The ciphertext.
 */
int glwe_secret_demasking(GLWECtParams* enc_params,
                          double* phase,
                          GGSWPreparedSK* sk_dft,
                          GLWECiphertext* ct);

/**
 * @brief Encrypts the phase (message + noise) and puts it in res.
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
int glwe_secret_masking(const MODULE* module,
                        GLWECtParams* params,
                        VecBiv* res_ct,
                        GGSWPreparedSK* sk,
                        PolyBiv* phase);


/**
 * @brief Encrypts the message m into GGSW ciphertext res with parameters enc_params.
 * 
 * @param enc_params The encryption params
 * @param res The encrypted message
 * @param sk_dft The secret key
 * @param msg_univ The message
 * 
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int ggsw_secret_encrypt(GGSWCtParams* enc_params,
                        GGSWCiphertext* res,           
                        GGSWPreparedSK* sk_dft,             
                        PolyUniv* msg_univ);

// TODO : the 8 functions below are not implemented

int ggsw_secret_masking(GGSWCiphertext* res,
                        GGSWPreparedSK* sk,
                        PolyUniv* msg_univ,     
                        GGSWCtParams* enc_params);

void ggsw_secret_demasking(double* res,   // result
                  GGSWPreparedSK* sk_dft, // secret key
                  GGSWCiphertext* ct      // ciphertext
);

// GGSWPublicKey is a struct encapsulating everything regarding the public
// key.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_public_encrypt(GGSWCiphertext* res,           // result
                         GGSWPublicKey* pk,             // secret key
                         int64_t* m,                // message
                         GGSWCtParams* enc_params  // parameters
);

// PartialGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_secret_encrypt(PartialGGSWCiphertext* res,       // result
                             GGSWSecretKey* sk,             // secret key
                             int64_t* m,                // message
                             GGSWCtParams* enc_params  // parameters
);

/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_public_encrypt(PartialGGSWCiphertext* res,       // result
                             GGSWPublicKey* pk,             // public key
                             int64_t* m,                // message
                             GGSWCtParams* enc_params  // parameters
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void halfggsw_decrypt(int64_t* res,       // result
                      GGSWSecretKey* sk,      // secret key
                      PartialGGSWCiphertext* ct  // ciphertext
);

/* Should it be in glwe.h since result is GLWE ? */
void ggsw_external_product(GLWECiphertext* res,  // result
                           GLWECiphertext* ct1,  // GLWE ciphertext
                           GGSWCiphertext* ct2   // GGSW ciphertext
);

void halfggsw_external_product(GLWECiphertext* res,     // result
                               GLWECiphertext* ct1,     // GLWE ciphertext
                               PartialGGSWCiphertext* ct2  // half GGSW ciphertext
);

//! GGSW IN DFT PART (begin)   

/**
 * @brief Adds a bivariate error to the bivariate phase, returns in DFT space.
 * 
 * @param enc_params The GLWE parameters.
 * @param res_dft The result bivariate phase in DFT space.
 * @param phase_dft The input phase.
 */
void add_error_dft(GLWECtParams* enc_params,
                  PolyBivDFT* res_dft,
                  PolyBivDFT* phase_dft);

/**
 * @brief Masks the phase (message + noise) in DFT space and puts it in res_ct.
 * 
 * @param enc_params The GLWE parameters.
 * @param module The module stocking the degree N.
 * @param res_dft The result ciphertext in DFT space. 
 * @param sk_dft The secret key in DFT space.
 * @param phase_dft message + error.
 * 
 * @retval `-1` if an error occurs.
 * @retval `0` otherwise.
 */
int glwe_secret_masking_dft(GLWECtParams* enc_params,
                         const MODULE* module,
                         VecBivDFT* res_dft,
                         GGSWPreparedSK* sk_dft,
                         PolyBivDFT* phase_dft);

/**
 * @brief Encrypts the message m into GGSW ciphertext res with parameters enc_params in DFT space.
 * 
 * @param enc_params The encryption params
 * @param res_dft The encrypted message
 * @param sk_dft The secret key
 * @param msg_univ The message
 * 
 * @retval `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval `0` otherwise.
 */
int ggsw_secret_encrypt_dft(GGSWCtParams* enc_params,
                            GGSWCiphertextDFT* res_dft,
                            GGSWPreparedSK* sk_dft,
                            PolyUniv* msg_univ);


//! COMMON PART (begin)

#endif  // GGSW_H
