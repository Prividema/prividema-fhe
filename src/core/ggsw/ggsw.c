#include "ggsw.h"
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
){

}

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
){
    // The number of terms
    int64_t k = ct1->params->k;

    for(int i=0; i <k; ++i){

        halfggsw_add(res->values + i , ct1->values + i , ct2->values + i);

    }
}

void ggsw_add_inplace(GGSWCiphertext* res,  // result
                      GGSWCiphertext* ct    // ciphertext
);

/* Should it be in glwe.h since result is GLWE ? */
void ggsw_external_product(GLWECiphertext* res,  // result
                           GLWECiphertext* ct1,  // GLWE ciphertext
                           GGSWCiphertext* ct2   // GGSW ciphertext
);

/** @brief Adds two halfggsw */
void halfggsw_add(HalfGGSWCiphertext* res,  // result
                  HalfGGSWCiphertext* ct1,  // Half GGSW ciphertext
                  HalfGGSWCiphertext* ct2  // Half GGSW ciphertext
){

    int64_t l_gadget = res->params->l;
    for(int64_t i = 0; i < l_gadget ; ++i){

    }
}

                  
void halfggsw_external_product(GLWECiphertext* res,     // result
                               GLWECiphertext* ct1,     // GLWE ciphertext
                               HalfGGSWCiphertext* ct2  // half GGSW ciphertext
                            
);