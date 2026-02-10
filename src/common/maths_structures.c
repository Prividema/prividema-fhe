#include "maths_structures.h"
#include <stdio.h>

void printf_poly_biv(PolyBiv* pol, int64_t pol_sl, int64_t N, int64_t l){
    printf("\n");
    for(int64_t i = 1 ; i <= l ; i++)
    {
        printf("\nY^%ld : ", i);
        for(int64_t p = 0 ; p < N ; p++){
            if(pol[(i-1)*pol_sl + p] < 0)
                printf("%ld X^%ld ", pol[(i-1)*pol_sl + p], p);
            else
                printf(" %ld X^%ld ", pol[(i-1)*pol_sl + p], p);
        }
    }
}

void printf_vec_poly_biv(VecBiv* pols, int64_t pols_nb, int64_t N, int64_t l){
    printf("\n");
    for(int64_t j = 0 ; j < pols_nb ; j++){
        printf("\n%ld-th component", j);
        printf_poly_biv(pols + j*N, pols_nb*N, N, l);
        printf("\n");
    }
}

void printf_poly_univ_ZnX(PolyUniv* pol, int64_t N){
    for(int64_t p = 0 ; p < N ; p++){
        printf(" %ld X^%ld ", pol[p], p);
    }
}

void printf_poly_univ_RnX(double* pol, int64_t N){
    for(int64_t p = 0 ; p < N ; p++){
        printf(" %lf X^%ld ", pol[p], p);
    }
}


void printf_vec_poly_univ(VecBiv* pols, int64_t pols_size, int64_t N){
    printf("\n");
    for(int64_t j = 0 ; j < pols_size ; j++){
        printf("\n%ld-th component", j);
        printf_poly_univ_ZnX(pols + j*N, N);
        printf("\n");
    }
}

void printf_secret_key(PolyUniv** sk_values, int64_t N, int64_t k){
    printf("\n\nBegin Secret Key:\n ");
    for(int64_t j = 0 ; j < k ; j++){
        printf("\n%ld-th component", j);
        printf_poly_univ_ZnX(sk_values[j], N);
    }
    printf("\n\nEnd Secret Key\n ");
}
