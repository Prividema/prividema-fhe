#ifndef GGSW_H
#define GGSW_H

#include <cstdint>
#include "ggsw_encrypt_params.h"
#include "ggsw_ciphertext.h"
#include "ggsw_secret_key.h"
#include "ggsw_public_key.h"
#include "halfggsw_ciphertext.h"
//#include "something that defines INT_POL and BASE2K_INT_POL.h"


// GSW_CIPHERTEXT is a struct encapsulating ciphertext values and params.
// GSW_SECRET_KEY is a struct encapsulating everything regarding the secret key.
// INT_POL is a struct encapsulating a polynomial with integer coefficients. Probably ZNX.
// GSW_ENCRYPT_PARAMS encapsulate all the encryption parameters.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void gsw_secret_encrypt(
    GSW_CIPHERTEXT* res,
    GSW_SECRET_KEY* sk, 
    INT_POL* m,
    GSW_ENCRYPT_PARAMS* enc_params
);

// GSW_PUBLIC_KEY is a struct encapsulating everything regarding the public key.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void gsw_public_encrypt(
    GSW_CIPHERTEXT* res,
    GSW_PUBLIC_KEY* pk,
    INT_POL* m,
    GSW_ENCRYPT_PARAMS* enc_params
);

// HALFGSW_CIPHERTEXT is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfgsw_secret_encrypt(
    HALFGSW_CIPHERTEXT* res,
    GSW_SECRET_KEY* sk,
    INT_POL* m,
    GSW_ENCRYPT_PARAMS* enc_params
);

/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfgsw_public_encrypt(
    HALFGSW_CIPHERTEXT* res,
    GSW_PUBLIC_KEY* pk,
    INT_POL* m,
    GSW_ENCRYPT_PARAMS* enc_params
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void gsw_decrypt(
    INT_POL* res,
    GSW_SECRET_KEY* sk,
    GSW_CIPHERTEXT* ct
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void halfgsw_decrypt(
    INT_POL* res,
    GSW_SECRET_KEY* sk,
    HALFGSW_CIPHERTEXT* ct
);

/* Adds two GSW ciphertext with same params and put result in res */
void gsw_add(
    GSW_CIPHERTEXT* res,
    GSW_CIPHERTEXT* ct1,
    GSW_CIPHERTEXT* ct2
);

/* Should it be in glwe.h since result is GLWE ? */
void gsw_external_product(
    GLWE_CIPHERTEXT* res,
    GLWE_CIPHERTEXT* ct1,
    GLWE_CIPHERTEXT* ct2
);

void halfGSW_external_product(
    GLWE_CIPHERTEXT* res,
    GLWE_CIPHERTEXT* ct1,
    HALFGSW_CIPHERTEXT* ct2
);

#endif
