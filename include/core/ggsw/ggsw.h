#ifndef GGSW_H
#define GGSW_H

#include "ggsw_encrypt_params.h"
#include "ggsw_key.h"
#include "glwe_ciphertext.h"
#include "halfggsw_ciphertext.h"

// #include "something that defines INT_POL and BASE2K_INT_POL.h"

// GGSWCiphertext is a struct encapsulating ciphertext values and params.
// GGSWSecretKey is a struct encapsulating everything regarding the secret
// key. INT_POL is a struct encapsulating a polynomial with integer
// coefficients. Probably ZNX. GGSWEncryptParams encapsulate all the
// encryption parameters.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_secret_encrypt(GGSWCiphertext* res,           // result
                         GGSWSecretKey* sk,             // secret key
                         IntegerPoly* m,                // message
                         GGSWEncryptParams* enc_params  // parameters
);

// GGSWPublicKey is a struct encapsulating everything regarding the public
// key.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_public_encrypt(GGSWCiphertext* res,           // result
                         GGSWPublicKey* pk,             // secret key
                         IntegerPoly* m,                // message
                         GGSWEncryptParams* enc_params  // parameters
);

// HalfGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_secret_encrypt(HalfGGSWCiphertext* res,       // result
                             GGSWSecretKey* sk,             // secret key
                             IntegerPoly* m,                // message
                             GGSWEncryptParams* enc_params  // parameters
);

/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_public_encrypt(HalfGGSWCiphertext* res,       // result
                             GGSWPublicKey* pk,             // public key
                             IntegerPoly* m,                // message
                             GGSWEncryptParams* enc_params  // parameters
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void ggsw_decrypt(IntegerPoly* res,   // result
                  GGSWSecretKey* sk,  // secret key
                  GGSWCiphertext* ct  // ciphertext
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void halfggsw_decrypt(IntegerPoly* res,       // result
                      GGSWSecretKey* sk,      // secret key
                      HalfGGSWCiphertext* ct  // ciphertext
);

/* Adds two GGSW ciphertext with same params and put result in res */
void ggsw_add(GGSWCiphertext* res,  // result
             GGSWCiphertext* ct1,  // first operand
             GGSWCiphertext* ct2   // second operand
);

void ggsw_add_inplace(GGSWCiphertext* res,  // result
                      GGSWCiphertext* ct    // ciphertext
);

/* Should it be in glwe.h since result is GLWE ? */
void ggsw_external_product(GLWECiphertext* res,  // result
                           GLWECiphertext* ct1,  // GLWE ciphertext
                           GGSWCiphertext* ct2   // GGSW ciphertext
);

void halfggsw_external_product(GLWECiphertext* res,     // result
                               GLWECiphertext* ct1,     // GLWE ciphertext
                               HalfGGSWCiphertext* ct2  // half GGSW ciphertext
);
#endif  // GGSW_H
