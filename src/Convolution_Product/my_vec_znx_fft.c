#include "my_vec_znx_fft.h"

void print_vec_znx(int64_t* p, int64_t L, int64_t N){
    for (int64_t i=0; i < L; i++){
        for (int64_t j=0; j < N; j++){
            printf(" %liX^%li ",p[ i*N + j ],j);
        }
        printf("\n");
    }
    printf("\n");
}
void print_vec_znx_dft(double* p, int64_t L, int64_t N){
    for (int64_t i=0; i < L; i++){
        for (int64_t j=0; j < N; j++){
            printf(" %f + i%f (%li) ",p[ 2*i*N + j ],p[ (2*i+1)*N + j], j);
        }
        printf("\n");
    }
    printf("\n");
}

/** @brief sets res = a+b */
void my_vec_dft_add(const MODULE* module,                   // N
                        VEC_ZNX_DFT* res, uint64_t res_size,    // res
                        const VEC_ZNX_DFT* a, uint64_t a_size,  // a
                        const VEC_ZNX_DFT* b, uint64_t b_size   // b
){

  const uint64_t m = module->m;
  if (a_size <= b_size) {
    const uint64_t sum_idx = res_size < a_size ? res_size : a_size;
    const uint64_t copy_idx = res_size < b_size ? res_size : b_size;
    // add up to the smallest dimension
    for (uint64_t i = 0; i < sum_idx; ++i) {

      for(uint64_t j = 0; j < m; ++j){

		// Re Part
        ((double *)res)[2*i*m + j] = ((double *)a)[2*i*m + j] + ((double *)b)[2*i*m + j];
		// Im Part
		((double *)res)[(2*i + 1)*m + j] = ((double *)a)[(2*i + 1)*m + j] + ((double *)b)[(2*i + 1)*m + j];

      }
    }
    // then copy to the largest dimension
    for (uint64_t i = sum_idx; i < copy_idx; ++i) {

        for(uint64_t j = 0; j<m; ++j){

            // Re Part
			((double *)res)[2*i*m + j] = ((double *)b)[2*i*m + j];
			// Im Part
			((double *)res)[(2*i + 1)*m + j] = ((double *)b)[(2*i + 1)*m + j];

      }

    }    
    // then extend with zeros
    for (uint64_t i = copy_idx; i < res_size; ++i) {
		
		uint64_t nn = module->nn;
		memset(((double*)res) + i*nn , 0 , nn*sizeof(double));
        
    }
  } else {

    my_vec_dft_add(module, res, res_size, b, b_size, a, a_size);

  }
}

/** @brief sets res = a-b */
void my_vec_dft_sub(const MODULE* module,                   // N
                        VEC_ZNX_DFT* res, uint64_t res_size,    // res
                        const VEC_ZNX_DFT* a, uint64_t a_size,  // a
                        const VEC_ZNX_DFT* b, uint64_t b_size   // b
){

  const uint64_t m = module->m;
  if (a_size <= b_size) {
    const uint64_t sum_idx = res_size < a_size ? res_size : a_size;
    const uint64_t copy_idx = res_size < b_size ? res_size : b_size;
    // sub up to the smallest dimension
    for (uint64_t i = 0; i < sum_idx; ++i) {

      for(uint64_t j = 0; j < m; ++j){

		// Re Part
        ((double *)res)[2*i*m + j] = ((double *)a)[2*i*m + j] - ((double *)b)[2*i*m + j];
		// Im Part
		((double *)res)[(2*i + 1)*m + j] = ((double *)a)[(2*i + 1)*m + j] - ((double *)b)[(2*i + 1)*m + j];

      }
    }
    // then copy to the largest dimension
    for (uint64_t i = sum_idx; i < copy_idx; ++i) {

        for(uint64_t j = 0; j<m; ++j){

            // Re Part
			((double *)res)[2*i*m + j] = -((double *)b)[2*i*m + j];
			// Im Part
			((double *)res)[(2*i + 1)*m + j] = -((double *)b)[(2*i + 1)*m + j];

      }

    }    
    // then extend with zeros
    for (uint64_t i = copy_idx; i < res_size; ++i) {
		
		uint64_t nn = module->nn;
		memset(((double*)res) + i*nn , 0 , nn*sizeof(double));
        
    }
  } else {
		my_vec_dft_sub(module, res, res_size, b, b_size, a, a_size);
  }
}

/** @brief In the k-eme term of the bivariate product a * b, Computes FFT(a_i*b_(k-i)) 
*/
void my_vec_dft_conv_k_i(const MODULE* module, 
                    VEC_ZNX_DFT* res, 
					const VEC_ZNX_DFT* a, const VEC_ZNX_DFT* b,
					int64_t k, int64_t i
){
	uint64_t m = module->m;
    /*
    printf("hi_ki %li %li %li\n", m,k, i);
    printf("a_fft : ");
	print_vec_znx_dft((double *)a, 1, m);
    printf("b_fft : ");
	print_vec_znx_dft((double *)b, 1, m);*/
	// Element-wise Product FFT(a_i*b_(k-i)) = FFt(a_i) * FFT(b_(i-k))
	for(int64_t p=0; p < m; ++p){
    double a_re = ((double*)a)[2*i*m + p];
	double a_im = ((double*)a)[(2*i + 1)*m + p];

	double b_re = ((double*)b)[2*(k-i)*m + p];
	double b_im = ((double*)b)[(2*(k-i) + 1)*m + p];
	
	// Real part of Res
	((double*)res)[2*k*m + p] += a_re * b_re - a_im * b_im;

	// Imaginary part of Res
	((double *)res)[(2*k + 1)*m + p] += a_re * b_im + a_im * b_re;
	}
}
/** @brief computes the k-eme term of the bivariate product FFT(a * b) */
void my_vec_dft_conv_k(const MODULE* module, 
    VEC_ZNX_DFT* res,
    const VEC_ZNX_DFT* a, const VEC_ZNX_DFT* b,
    const int64_t k, uint64_t min_i, uint64_t max_i
){
    for (int64_t i = min_i ; i <= max_i; ++i){
        printf("hi_k %li %li %li\n", i, min_i, max_i);
	    my_vec_dft_conv_k_i(module, res, a, b, k, i);

	}
}
  

/** @brief Computes the convolution product of two bivariate polynomial in DFT Space */
void my_vec_dft_conv(const MODULE* module, 
                    VEC_ZNX_DFT* res, uint64_t res_size,
                    const VEC_ZNX_DFT* a, uint64_t a_size,
                    const VEC_ZNX_DFT* b, uint64_t b_size
){	
	uint64_t conv_size = a_size + b_size - 1;

	if (a_size <= b_size){

		uint64_t conv_fst_idx = a_size < res_size ? a_size : res_size;
		uint64_t conv_snd_idx = b_size < res_size ? b_size : res_size;
		uint64_t conv_thd_idx = conv_size < res_size ? conv_size : res_size;
		printf("hi_conv: %li %li %li\n",conv_fst_idx,conv_snd_idx,conv_thd_idx);
		for(int64_t k = 0 ; k < conv_fst_idx ; ++k){
			my_vec_dft_conv_k(module, res, a, b, k, 0, k);
		}
		for(int64_t k = conv_fst_idx ; k < conv_snd_idx ; ++k){
            printf("holla");
			my_vec_dft_conv_k(module, res, a, b, k, 0, a_size);
		}
		for(int64_t k = conv_snd_idx ; k < conv_thd_idx ; ++k){
            printf("holla");
			my_vec_dft_conv_k(module, res, a, b, k, k-b_size, a_size);
		}
		for(int64_t i = conv_thd_idx ; i <res_size ; ++i){
            printf("holla");
			uint64_t nn = module->nn;
			memset(((double*)res) + i*nn , 0 , nn*sizeof(double));
		}
	}

	else {
		my_vec_dft_conv(module, res, res_size, b, b_size, a, a_size);
	}
}
