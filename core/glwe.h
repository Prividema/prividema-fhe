#include <cassert>
#include <cstdint>
#include <iostream>
#include <cstdint>
#include <random>

struct glwe_ciphertext {
    uint64_t N; 
    uint64_t k;
    uint64_t kappa; //base_2k
    uint64_t l_a; //size of the decomposition of a
    uint64_t l_b; //size of the decomposition of b
    int64_t* b; 
    int64_t* a; // vector of k elements
    void* data; // memory of the struct
}

struct glwe_key {
    uint64_t n; // n= N*k (n is constant: N and k can change when we change the module N) 
    int64_t* s; 
    void* data; 
}

struct glwe_prep_key {
    uint64_t N;
    uint64_t k;
    SVP_PPOL** s; // vec of size k, each element is prepared vec
    void* data;
}

//tnx, rthild....
struct tnx_element {
    uint64_t N;
    uint64_t kappa;
    uint64_t l;
    int64_t* coeffs;
    void* data;
}

//encrypt
void glwe_encrypt_fft(  
                         const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used   
                         glwe_ciphertext* ct,                                      //ciphertext
                         glwe_prep_key* s,                                         // secret key: vec of size k
                         tnx_element* phase                                        // message + noise
);



//decrypt (compute the phase)
void glwe_phase(       const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used  
                       tnx_element* phase,                                      // decrypted phase
                       glwe_prep_key* s,                                        // secret key
                       glwe_ciphertext* ct,                                      //ciphertext
);




