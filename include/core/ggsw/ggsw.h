#ifndef GGSW_H
#define GGSW_H

#include "ggsw_encrypt_params.h"
#include "ggsw_key.h"
#include "halfggsw_ciphertext.h"
#include "glwe_ciphertext.h"

//#include "something that defines INT_POL and BASE2K_INT_POL.h"

// GGSW_CIPHERTEXT is a struct encapsulating ciphertext values and params.
// GGSW_SECRET_KEY is a struct encapsulating everything regarding the secret key.
// INT_POL is a struct encapsulating a polynomial with integer coefficients. Probably ZNX.
// GGSW_ENCRYPT_PARAMS encapsulate all the encryption parameters.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_secret_encrypt(
    GGSW_CIPHERTEXT* res,
    GGSW_SECRET_KEY* sk, 
    INT_POL* m,
    GGSW_ENCRYPT_PARAMS* enc_params
);

// GGSW_PUBLIC_KEY is a struct encapsulating everything regarding the public key.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_public_encrypt(
    GGSW_CIPHERTEXT* res,
    GGSW_PUBLIC_KEY* pk,
    INT_POL* m,
    GGSW_ENCRYPT_PARAMS* enc_params
);

// HALFGSW_CIPHERTEXT is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_secret_encrypt(
    HALFGGSW_CIPHERTEXT* res,
    GGSW_SECRET_KEY* sk,
    INT_POL* m,
    GGSW_ENCRYPT_PARAMS* enc_params
);

/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_public_encrypt(
    HALFGGSW_CIPHERTEXT* res,
    GGSW_PUBLIC_KEY* pk,
    INT_POL* m,
    GGSW_ENCRYPT_PARAMS* enc_params
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void ggsw_decrypt(
    INT_POL* res,
    GGSW_SECRET_KEY* sk,
    GGSW_CIPHERTEXT* ct
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void halfggsw_decrypt(
    INT_POL* res,
    GGSW_SECRET_KEY* sk,
    HALFGGSW_CIPHERTEXT* ct
);


/* Adds two GGSW ciphertext with same params and put result in res */
void gsw_add(
    GGSW_CIPHERTEXT* res,
    GGSW_CIPHERTEXT* ct1,
    GGSW_CIPHERTEXT* ct2
);

void ggsw_add_inplace(
    GGSW_CIPHERTEXT* res,
    GGSW_CIPHERTEXT* ct
);

/* Should it be in glwe.h since result is GLWE ? */
void ggsw_external_product(
    GLWE_CIPHERTEXT* res,
    GLWE_CIPHERTEXT* ct1,
    GGSW_CIPHERTEXT* ct2
);

void halfggsw_external_product(
    GLWE_CIPHERTEXT* res,
    GLWE_CIPHERTEXT* ct1,
    HALFGGSW_CIPHERTEXT* ct2
);
#endif // GGSW_H