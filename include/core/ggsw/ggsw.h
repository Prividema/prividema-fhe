#ifndef GGSW_H
#define GGSW_H

#include "ggsw_ct_params.h"
#include "ggsw_key.h"
#include "glwe_ciphertext.h"
#include "partialggsw_ciphertext.h"

#define WITH_Y0

void add_vec_znxy(GGSWCiphertext* res, GGSWCiphertext* c, GGSWCiphertext* d);

// GGSWCiphertext is a struct encapsulating ciphertext values and params.
// GGSWSecretKey is a struct encapsulating everything regarding the secret
// key. INT_POL is a struct encapsulating a polynomial with integer
// coefficients. Probably ZNX. GGSWCtParams encapsulate all the
// encryption parameters.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_secret_encrypt(GGSWCiphertext* res,           // result
                         GGSWSecretKey* sk,             // secret key
                         int64_t* m,                // message
                         GGSWCtParams* enc_params  // parameters
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

/* Adds two GGSW ciphertext with same params and put result in res */
void gsw_add(GGSWCiphertext* res,  // result
             GGSWCiphertext* ct1,  // first operand
             GGSWCiphertext* ct2   // second operand
);

int add(int a, int b);

int multiply(int a, int b);

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
                               PartialGGSWCiphertext* ct2  // half GGSW ciphertext
);
#endif  // GGSW_H
