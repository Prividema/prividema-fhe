#include "ggsw.h"

/** @brief Print c_ggsw in FT space */
void print_ggsw_fft(Core* core, GGSWCiphertext* c_ggsw){
    uint64_t k = core->k;
    uint64_t N = core->R->nn;

}

/** @brief Print c_ggsw in ZNX space */
void print_ggsw(Core* core, GGSWCiphertext* c_ggsw){
    uint64_t k = core->k;
    uint64_t N = core->R->nn;

}

/** @brief Print c_halfggsw in FT space */
void print_halfggsw_fft(Core* core, GGSWCiphertext* c_halfggsw){
    uint64_t k = core->k;
    uint64_t N = core->R->nn;

}

/** @brief Print c_halfggsw in ZNX space */
void print_halfggsw(Core* core, HalfGGSWCiphertext* c_halfggsw){
    uint64_t k = core->k;
    uint64_t N = core->R->nn;
    uint64_t degreeInY = c_halfggsw->params->DegreeInY;
    uint64_t nrows = c_halfggsw->params->nrows;

    for(int i = 0 ; i < nrows ; i++){

        for(int j = 0 ; j < (k+1) ; j++){

        }
    }
}
// GGSWCiphertext is a struct encapsulating ciphertext values and params.
// GGSWSecretKey is a struct encapsulating everything regarding the secret
// key. INT_POL is a struct encapsulating a polynomial with integer
// coefficients. Probably ZNX. GGSWEncryptParams encapsulate all the
// encryption parameters.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_secret_encrypt(Core* core,
                         GGSWCiphertext* res,           // result
                         GGSWSecretKey* sk,             // secret key
                         IntegerPoly* m,                // message
                         GGSWEncryptParams* enc_params  // parameters
){

}

// GGSWPublicKey is a struct encapsulating everything regarding the public
// key.
/* Encrypts message m into GGSW ciphertext res with parameters enc_params */
void ggsw_public_encrypt(Core* core, 
                         GGSWCiphertext* res,           // result
                         GGSWPublicKey* pk,             // secret key
                         IntegerPoly* m,                // message
                         GGSWEncryptParams* enc_params  // parameters
){

}

// HalfGGSWCiphertext is a struct encapsulating ciphertext values and params.
/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_secret_encrypt(Core* core, HalfGGSWCiphertext* res,       // result
                             GGSWSecretKey* sk,             // secret key
                             IntegerPoly* m,                // message
                             GGSWEncryptParams* enc_params  // parameters
);

/* Encrypts message m into halfGGSW ciphertext res with parameters enc_params */
void halfggsw_public_encrypt(Core* core, HalfGGSWCiphertext* res,       // result
                             GGSWPublicKey* pk,             // public key
                             IntegerPoly* m,                // message
                             GGSWEncryptParams* enc_params  // parameters
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void ggsw_decrypt(Core* core, IntegerPoly* res,   // result
                  GGSWSecretKey* sk,  // secret key
                  GGSWCiphertext* ct  // ciphertext
);

/* Decrypts message res from GGSW ciphertext ct */
/* /!\ Is that actually useful ? */
void halfggsw_decrypt(Core* core, IntegerPoly* res,       // result
                      GGSWSecretKey* sk,      // secret key
                      HalfGGSWCiphertext* ct  // ciphertext
);


/* Should it be in glwe.h since result is GLWE ? */
void ggsw_external_product(Core* core, GLWECiphertext* res,  // result
    GLWECiphertext* ct1,  // GLWE ciphertext
    GGSWCiphertext* ct2   // GGSW ciphertext
);

/** @brief Adds two halfggsw */
void halfggsw_add(Core* core, HalfGGSWCiphertext* res,  // result
    HalfGGSWCiphertext* ct1,  // Half GGSW ciphertext
    HalfGGSWCiphertext* ct2  // Half GGSW ciphertext
){
    /** FOR NOW UNABLE TO ADD IF DIFFERENT SIZE ? */
    if ((ct1->params->nrows != ct2->params->nrows) || (ct1->params->nrows != ct2->params->nrows)){
        printf("COMPILATION ERROR : CT1 and CT2 are not the same size !");
    }
    
    uint64_t nrows = res->params->nrows;
    uint64_t k = core->k;
    uint64_t N = core->R->nn;
    uint64_t MinDegreeInY = ct1->params->DegreeInY < ct2->params->DegreeInY ? 
                          ct1->params->DegreeInY : ct2->params->DegreeInY ;

    uint64_t* res_mat = (uint64_t*)res->values;
    uint64_t* ct1_mat = (uint64_t*)ct1->values;
    uint64_t* ct2_mat = (uint64_t*)ct2->values;

    for(int64_t i = 0 ; i < nrows ; ++i){
        
        for(int64_t j = 0 ; j < k+1 ; ++j){

            for(int64_t jj = 0 ; jj < MinDegreeInY ; ++jj){

                for(int64_t jjj = 0 ; jjj < N ; ++jjj){

                    uint64_t index = res->index(k, MinDegreeInY, N, i, j, jj, jjj);
                    res_mat[index] = ct1_mat[index] + ct2_mat[index];

                }
            }
        }
    }
}

/* Adds two GGSW ciphertext with same params and put result in res */
void ggsw_add(Core* core, GGSWCiphertext* res,  // result
             GGSWCiphertext* ct1,  // first operand
             GGSWCiphertext* ct2   // second operand
){
    // The number of terms
    int64_t k = core->k;

    for(int i=0; i <k+1; ++i){

        halfggsw_add(core, res->halfggswVec + i , ct1->halfggswVec + i , ct2->halfggswVec + i);

    }
}

void ggsw_add_inplace(Core* core, GGSWCiphertext* res,  // result
                      GGSWCiphertext* ct    // ciphertext
);

void halfggsw_external_product(Core* core, GLWECiphertext* res,     // result
    GLWECiphertext* ct1,     // GLWE ciphertext
    HalfGGSWCiphertext* ct2  // half GGSW ciphertext
                            
);