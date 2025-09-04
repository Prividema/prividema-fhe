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
struct rthild {
    uint64_t N;
    uint64_t kappa;
    uint64_t l;
    int64_t* coeffs;
    void* data;
}


void glwe_encrypt_fft(  
                         const CORE* core,                                        // all params of the library: is fft or ntt, all N that are used   
                         uint64_t log2_base2k,                                    // output base 2^Kappa (Kappa is number of bits in the decomposition for each element)
                         uint64_t k;                                               // number of samples (k=1 when RLWE)
                         int64_t* new_b, uint64_t new_b_size, uint64_t new_b_sl,  // b part of glwe
                         const SVP_PPOL* s,                                       // secret key: vec of size k
                         const int64_t* a, uint64_t a_size, uint64_t a_sl,        // a part of glwe
                         const int64_t* phi, uint64_t phi_size, uint64_t phi_sl   // message + noise
);




void rlwe_phase(const MODULE* module,                              //
                       uint64_t log2_base2k,                              // output base 2^K
                       int64_t* phi, uint64_t phi_size, uint64_t phi_sl,  // decrypted phase
                       const SVP_PPOL* s,                                 // secret key
                       const int64_t* a, uint64_t a_size, uint64_t a_sl,  // a part of rlwe
                       const int64_t* b, uint64_t b_size, uint64_t b_sl   // message + noise
);




