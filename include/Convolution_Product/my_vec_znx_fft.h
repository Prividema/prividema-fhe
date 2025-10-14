#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "vec_znx_arithmetic_private.h" 

/** @brief Print a vec ZnX */
void print_vec_znx(int64_t* p, int64_t L, int64_t N);

/** @brief Print a vec ZnX in DFT space*/
void print_vec_znx_dft(double* p, int64_t L, int64_t N);

/**  */
void my_vec_dft_add(const MODULE* module,                   // N
                    VEC_ZNX_DFT* res, uint64_t res_size,    // res
                    const VEC_ZNX_DFT* a, uint64_t a_size,  // a
                    const VEC_ZNX_DFT* b, uint64_t b_size   // b
);
void my_vec_dft_sub(const MODULE* module, 
                    VEC_ZNX_DFT* res, uint64_t res_size,
                    const VEC_ZNX_DFT* a, uint64_t a_size,
                    const VEC_ZNX_DFT* b, uint64_t b_size
);
void my_vec_dft_conv_k_i(const MODULE* module, 
                    VEC_ZNX_DFT* res, 
					const VEC_ZNX_DFT* a, const VEC_ZNX_DFT* b,
					int64_t k, int64_t i
);
void my_vec_dft_conv_k(const MODULE* module, 
                    VEC_ZNX_DFT* res,
                    const VEC_ZNX_DFT* a, const VEC_ZNX_DFT* b,
                    const int64_t k, uint64_t min_i, uint64_t max_i
);
void my_vec_dft_conv(const MODULE* module, 
                    VEC_ZNX_DFT* res, uint64_t res_size,
                    const VEC_ZNX_DFT* a, uint64_t a_size,
                    const VEC_ZNX_DFT* b, uint64_t b_size);
