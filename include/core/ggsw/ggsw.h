#ifndef GGSW_H
#define GGSW_H

#include "ggsw_ciphertext.h"
#include "ggsw_key.h"
#include "spqlios_alias.h"

#define WITH_Y0


//! GGSW PART (begin)

int encrypt_biv_glwe(const MODULE* module,
                     GLWECtParams* params, 
                     VecBiv* res_ct,
                     GGSWPreparedSK* sk, 
                     PolyBiv* phase
);

int compute_phase_biv(GGSWCtParams* params, 
                      PolyBiv* res, 
                      PolyUniv* phase_univ, 
                      int64_t i_tilde
);

int ggsw_secret_encrypt(GGSWCiphertext* res,           
                        GGSWPreparedSK* sk,             
                        PolyUniv* msg_univ,                
                        GGSWCtParams* enc_params 
);

void ggsw_decrypt(double* res,   // result
                  GGSWPreparedSK* sk_dft,  // secret key
                  GGSWCiphertext* ct  // ciphertext
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

int decrypt_biv_glwe(GLWECtParams* enc_params,
                     double* phase, 
                     GGSWPreparedSK* sk_dft,
                     GLWECiphertext* ct
);

//! GGSW IN DFT PART (begin)   

int compute_phase_biv_dft(GGSWCtParams* enc_params, 
                          PolyBivDFT* res_dft,
                          PolyUnivDFT* phase_univ_dft, 
                          int64_t i
);

int add_error_dft(GLWECtParams* enc_params,
                  PolyBivDFT* res_dft,
                  PolyBivDFT* phase_dft 
);

int encrypt_biv_glwe_dft(GLWECtParams* enc_params, 
                         const MODULE* module, 
                         VecBivDFT* res_dft,
                         GGSWPreparedSK* sk_dft, 
                         PolyBivDFT* phase_dft
);

int ggsw_secret_encrypt_dft(GGSWCtParams* enc_params,    // parameters
                            GGSWCiphertextDFT* res_dft,        // result
                            GGSWPreparedSK* sk_dft,         // secret key
                            PolyUniv* msg_univ                 // message                
);


//! COMMON PART (begin)

int* add(int* a, int a_size, int* b, int b_size);
int add_int(int a, int b);
int multiply(int a, int b);
#endif  // GGSW_H
