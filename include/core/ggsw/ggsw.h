#ifndef GGSW_H
#define GGSW_H

#include "ggsw_ct_params.h"
#include "ggsw_key.h"
#include "glwe_ciphertext.h"
#include "partialggsw_ciphertext.h"

#define WITH_Y0

/** @brief opaque type that are included in */
typedef struct module_info_t MODULE;
/** @brief opaque type that represents a vector of znx in DFT space */
typedef struct vec_znx_dft_t VEC_ZNX_DFT;

void dft_mult_add_j0k_i_inplace(GLWECtParams* params,
                        double* a_j0k_i_dft,  
                        double* sk_dft
);



void secret_key_mult_add_inplace(GLWECtParams* params,
                         double* res_ct_dft,
                         double* sk_dft
);

void add_inplace_m_to_ap_dft(GLWECtParams* params, 
                         double* ct_in,
                         double* m,
                         int64_t p
);

int encrypt_biv_glwe(GLWECtParams* params, 
                     int64_t* res_ct,
                     const MODULE* module, 
                     VEC_ZNX_DFT* sk_dft, int64_t sk_size, 
                     VEC_ZNX_DFT* phase
);

void vec_znx_dft_mult(const MODULE* module, 
              VEC_ZNX_DFT* res_dft, int64_t res_size,
              VEC_ZNX_DFT* c_dft, int64_t c_size,  
              VEC_ZNX_DFT* d_dft, int64_t d_size
);


// GGSWCiphertext is a struct encapsulating ciphertext values and params.
// GGSWSecretKey is a struct encapsulating everything regarding the secret
// key. INT_POL is a struct encapsulating a polynomial with integer
// coefficients. Probably ZNX. GGSWCtParams encapsulate all the
// encryption parameters.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
int ggsw_secret_encrypt(GGSWCiphertext* res,           // result
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
                               PartialGGSWCiphertext* ct2  // half GGSW ciphertext
);

// Does the multiplication of (a_i)'s by the secret key (s_i)'s 
// We have : 
// 
// Result = ∑_i{0,k-1} s_i ⋅ a_i 
// Result = ∑_j{1,l} [ ∑_i{0,k-1} s_i ⋅ a_(i,j) ] ⋅ Y^j
void biv_secret_key_mult(int64_t N, int64_t k, int64_t l,
                         int64_t* res, 
                         int64_t* c,
                         int64_t* s
);

void decrypt_biv_glwe(int64_t* res, 
                      GGSWSecretKey* key,
                      int64_t* phase);

void ggsw_addition(GGSWCiphertext* ct_out,
                   GGSWCiphertext* ct_in1, GGSWCiphertext* ct_in2);


int* add(int* a, int a_size, int* b, int b_size);
int add_int(int a, int b);
int multiply(int a, int b);
#endif  // GGSW_H
