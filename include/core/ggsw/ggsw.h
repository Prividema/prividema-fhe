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

int ggsw_secret_encrypt(GGSWCiphertext* res,           
                        GGSWPreparedSK* sk,             
                        PolyUniv* msg,                
                        GGSWCtParams* enc_params 
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
void ggsw_decrypt(int64_t* res,   // result
                  GGSWSecretKey* sk,  // secret key
                  GGSWCiphertext* ct  // ciphertext
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

void decrypt_biv_glwe(int64_t* res, 
                      GGSWSecretKey* key,
                      int64_t* phase);


//! GGSW IN DFT PART (begin)   

int encrypt_biv_glwe_dft(GLWECtParams* params, 
                         VecBivDFT* res_ct,
                         const MODULE* module, 
                         GGSWPreparedSK* sk, 
                         PolyBivDFT* phase
);

int ggsw_secret_encrypt_dft(GGSWCiphertextDFT* res,        // result
                            GGSWPreparedSK* sk,         // secret key
                            PolyUniv* msg,                 // message
                            GGSWCtParams* enc_params    // parameters
);


//! COMMON PART (begin)

int* add(int* a, int a_size, int* b, int b_size);
int add_int(int a, int b);
int multiply(int a, int b);
#endif  // GGSW_H 
