#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "vec_znx_arithmetic_private.h" 
#include "my_vec_znx_fft.h"
int main(void) {
    // Parameters 
    int64_t N = 8;
    int64_t L = 2;

    // Module 
    MODULE* N_M = new_module_info(N,FFT64);

    // Bivariate polys
    int64_t* a = malloc(N*L*sizeof(int64_t));
    int64_t* b = malloc(N*L*sizeof(int64_t));
    VEC_ZNX_BIG* dest = new_vec_znx_big(N_M,2*L-1);

    // Init of a and b
    a[0] = 1;
    a[1] = 3;
	a[1 + 8] = 8;
	a[2 + 8] = 1;
    //a[3] = 4;
    b[0] = 1;
    b[5] = 1;
	b[1 + 8] = 8;
	b[2 + 8] = 1;
    //b[3] = 5;
	printf("a : ");
    print_vec_znx(a, L, N);
	printf("b : ");
    print_vec_znx(b, L, N);

    // FFT of those polys
    VEC_ZNX_DFT* a_fft = new_vec_znx_dft(N_M,L);
    VEC_ZNX_DFT* b_fft = new_vec_znx_dft(N_M,L);
    VEC_ZNX_DFT* dest_fft = new_vec_znx_dft(N_M,2*L-1);
    vec_znx_dft(N_M, a_fft, L, a, L, N);
    vec_znx_dft(N_M, b_fft, L, b, L, N);

    printf("a_fft : ");
	print_vec_znx_dft((double *)a_fft, L, N/2);
    printf("b_fft : ");
	print_vec_znx_dft((double *)b_fft, L, N/2);
    // FFT ADD
    //my_vec_dft_add(N_M,dest_fft, L,a_fft,L,b_fft,L);
    // FFT SUB

    // FFT CONV
	/*
	for (int i =0; i < N/2 ; ++i){
		double a = ((double *)a_fft)[i];
		double b = ((double *)a_fft)[i + N/2];
		double c = ((double *)b_fft)[i];
		double d = ((double *)b_fft)[i + N/2];
		((double*)dest_fft)[i] = a * c - b * d; 
		((double*)dest_fft)[i + N/2] = c * b + a * d;
	}*/
	my_vec_dft_conv(N_M, dest_fft, 2*L-1, a_fft, L, b_fft, L);
    // Result in DFT Space
    printf("dest_fft: ");
	print_vec_znx_dft(dest_fft,2*L-1,N/2);

    // Result in ZNX Space

    uint8_t* tmp;
    vec_znx_idft(N_M, dest, 2*L-1, dest_fft, 2*L-1, tmp );
    print_vec_znx(dest, 2*L-1, N);

    return 0;
}
