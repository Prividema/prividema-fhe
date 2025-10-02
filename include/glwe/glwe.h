#ifndef GLWE_H
#define GLWE_H

#include "glwe_key.h"

//tnx, rthild....
typedef struct tnx_element {
    uint64_t N;
    uint64_t kappa;
    uint64_t l;
    int64_t* coeffs;
    void* data;
} TNX_ELEMENT;

//secret key encrypt
void glwe_encrypt_priv(  
                         const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used   
                         GLWE_CIPHERTEXT* ct,                                     //ciphertext
                         GLWE_PREP_SECRET_KEY* s,                                        // secret key: vec of size k
                         TNX_ELEMENT* phase                                       // message + noise
);

//secret key decrypt (compute the phase)
void glwe_phase_priv(       const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used  
                       TNX_ELEMENT* phase,                                      // decrypted phase
                       GLWE_PREP_SECRET_KEY* s,                                        // secret key
                       GLWE_CIPHERTEXT* ct                                        //ciphertext
);

//public key encrypt
void glwe_encrypt_pub(  
                         const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used   
                         GLWE_CIPHERTEXT* ct,                                     //ciphertext
                         GLWE_PUBLIC_KEY* pk,                                     // public key  
                         TNX_ELEMENT* phase                                       // message + noise
);

//public key decrypt
void glwe_phase_pub(  
                         const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used   
                         GLWE_CIPHERTEXT* ct,                                     //ciphertext
                         GLWE_PREP_SECRET_KEY* s,                                        // secret key: vec of size k
                         TNX_ELEMENT* phase                                       // message + noise
);

//addition 2 glwe
void glwe_addition (
                       const CORE* core,  
                       GLWE_CIPHERTEXT* ct_out, 
                       GLWE_CIPHERTEXT* ct_in1,
                       GLWE_CIPHERTEXT* ct_in2
);

#endif // GLWE_H