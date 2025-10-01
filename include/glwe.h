#include <assert.h>
#include <stdint.h>

// TODO : To be Defined.
typedef struct core {} CORE;
typedef struct svp_ppol {} SVP_PPOL;

typedef struct glwe_ciphertext {
    uint64_t N;
    uint64_t k;
    uint64_t kappa; //base_2k
    uint64_t l_a; //size of the decomposition of a
    uint64_t l_b; //size of the decomposition of b
    int64_t* b; 
    int64_t* a; // vector of k elements
    void* data; // memory of the struct
} GLWE_CIPHERTEXT;

typedef struct glwe_secret_key {
    uint64_t n; // n= N*k (n is constant: N and k can change when we change the module N) 
    int64_t* s; 
    void* data; 
} GLWE_SECRET_KEY;

typedef struct glwe_prep_secret_key {
    uint64_t N;
    uint64_t k;
    SVP_PPOL** s; // vec of size k, each element is prepared vec
    void* data;
} GLWE_PREP_SECRET_KEY;

typedef struct glwe_public_key {
    uint64_t N;
    uint64_t k;
    uint64_t Y;
    GLWE_CIPHERTEXT** pk; // vector of Y element (A1,...,Ak, B)
    void* data; 
} GLWE_PUBLIC_KEY;

//tnx, rthild....
typedef struct tnx_element {
    uint64_t N;
    uint64_t kappa;
    uint64_t l;
    int64_t* coeffs;
    void* data;
} TNX_ELEMENT;

//generation secret key
void glwe_secret_key_gen (
                        const CORE* core,
                        uint64_t lambda, 
                        GLWE_SECRET_KEY* s
);

//generation public key
void glwe_public_key_gen (
                        const CORE* core,
                        uint64_t lambda, 
                        GLWE_SECRET_KEY* s,
                        GLWE_PUBLIC_KEY* pk
);

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
