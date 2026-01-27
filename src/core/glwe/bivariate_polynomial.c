#include "bivariate_polynomial.h"
#include "distributions.h"
#include "vec_znx_arithmetic_private.h"
#include <math.h>
#include <stdio.h>
#include <string.h>


//! BIV POLY PART (begin) 


uint64_t poly_biv_coef_number(GLWECtParams* params){
    uint64_t N = params->N;
    return poly_biv_size(params) * N;
}


PolyBiv* new_normal_random_biv_poly(MODULE* module, 
                                    GLWECtParams* params
){
    // GLWE parameters
    uint64_t N = params->N;
    uint64_t kappa = params->kappa;
    uint64_t l = params->n_limbs/(params->k + 1);

    // Draws a random univariate polynomial P(X) in Rn[X]
    double* rd_pol_univ = malloc(poly_univ_bytes(params));
    if(rd_pol_univ == NULL){
        perror("Malloc failed.");
        return NULL;
    }
    printf("%e\n", params->sigma);
    for(int64_t p = 0 ; p < N ; p++){
        if(rand_normal(rd_pol_univ + p, 0.0, params->sigma) < 0) {
            free(rd_pol_univ);
            return NULL;
        }
    }

    printf("rd_pol_univ : ");
    for(int64_t p = 0 ; p < N ; p++)
        printf("%e X^%ld ", rd_pol_univ[p], p);
    printf("\n");
    

    // Stores the base-2kappa normalized bivariate form Pbiv(X,Y) of P(X)
    PolyBiv* rd_pol = malloc(poly_biv_bytes(params));
    if(rd_pol == NULL){
        perror("Malloc failed.");
        free(rd_pol_univ);
        return NULL;
    }
    
    // For each, (p,i) in [0,N-1]x[0,l-1], Pbiv_p_i = centered(floor(P_p * Bg^i))
    // Where centered(_) is in [-2^(kappa-1) ; 2^(kappa-1) - 1]
    int64_t mask = (1LL << kappa) - 1;
    for(int64_t p = 0 ; p < N ; p++)
    {
        // For each p, we substract Bg/2 * Bg^(-i) to P_p
        for(int64_t i = 1 ; i < l ; i++){
            rd_pol_univ[p] += ldexp(1.0, kappa - 1 - kappa*i);
        }

        if(rd_pol_univ[p] >= 0)
        {
            for(int64_t i = 0 ; i < l ; i++)
            {
                rd_pol[i*N + p] = (((int64_t)ldexp(rd_pol_univ[p], i*kappa)) & mask) - (1LL << (kappa - 1));
            }
        }
        else{
            rd_pol_univ[p] -= floor(rd_pol_univ[p]);
            for(int64_t i = 0 ; i < l ; i++)
            {
                rd_pol[i*N + p] = (((int64_t)ldexp(rd_pol_univ[p], i*kappa)) & mask) - (1LL << (kappa - 1));
            }
        }
    }

    free(rd_pol_univ);

    return rd_pol;
}

void add_biv_poly(GLWECtParams* params, 
                  PolyBiv* res, int64_t res_sl,
                  PolyBiv* a, int64_t a_sl,
                  PolyBiv* b, int64_t b_sl 
){
    for(int64_t i = 0 ; i < poly_biv_size(params) ; i++)
    {
        for(int64_t p = 0 ; p < params->N ; p++)
        {
            res[p + i*res_sl] = a[p + i*a_sl] + b[p + i*b_sl];
        }
    }
}



//! BIV POLY IN DFT PART (begin) 


uint64_t poly_biv_coef_number_dft(GLWECtParams* params){
    uint64_t N = params->N;
    return (poly_biv_size(params) * N)/2;
}

PolyBivDFT* new_normal_random_biv_poly_dft(MODULE* module, 
                                           GLWECtParams* params
){
    // Base-2Kappa normalized bivariate polynomial in DFt space
    PolyBivDFT* rd_pol_dft = malloc(poly_biv_bytes(params));
    if(rd_pol_dft == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    // Base-2Kappa normalized bivariate polynomial
    PolyBiv* rd_pol = new_normal_random_biv_poly(module, params);
    if(rd_pol == NULL){
        perror("Malloc failed.");
        free(rd_pol_dft);
        return NULL;
    }
    
    // Then compute in DFT space
    vec_znx_dft_p(module, rd_pol_dft, poly_biv_size(params), rd_pol, poly_biv_size(params), params->N);

    free(rd_pol);
    
    return rd_pol_dft;
}

void add_biv_poly_dft(GLWECtParams* params, 
                      PolyBivDFT* res, int64_t res_sl,
                      PolyBivDFT* a, int64_t a_sl,
                      PolyBivDFT* b, int64_t b_sl 
){
    for(int64_t i = 0 ; i < poly_biv_size(params) ; i++)
    {
        for(int64_t p = 0 ; p < params->N ; p++)
        {
        res[p + i*res_sl] = a[p + i*a_sl] + b[p + i*b_sl];
        }
    }
    
}


//! COMMON PART (begin)

uint64_t poly_biv_bytes(GLWECtParams* params){
    return poly_biv_coef_number(params) * sizeof(int64_t);
}

uint64_t poly_biv_size(GLWECtParams* params){
    return params->n_limbs/(params->k + 1);
}

uint64_t poly_univ_bytes(GLWECtParams* params){
    uint64_t N = params->N;
    return N * sizeof(int64_t);
}

void biv_to_univ(GLWECtParams* params, double* pol_univ, PolyBiv* pol_biv){
    //GLWE parameters
    uint64_t N = params->N;
    uint64_t kappa = params->kappa;
    uint64_t l = poly_biv_size(params);

    // res_univ(X^p) = Sum_i{1,l}[poly(X^p, Y^i) * 2^(-kappa*i)]
    for(int64_t i = 1 ; i < l ; i++){
        for(int64_t p = 0 ; p < N ; p++){
            pol_univ[p] += ldexp((double)pol_biv[i*N + p], - i*kappa);
        }
    }
}

int univ_to_biv(GLWECtParams* params, PolyBiv* pol_biv, double* pol_univ
){
    // GLWE parameters
    uint64_t N = params->N;
    uint64_t kappa = params->kappa;
    uint64_t l = poly_biv_size(params);

    // Fills each pol_biv(X^p, Y^i) with the pol_univ's decomposition coefficients of  in [-2^(kappa* - 1) ; 2^(kappa - 1) - 1]
    int64_t mask = (1LL << kappa) - 1;

    double* tmp_pol_univ = calloc(poly_biv_bytes(params),1);

    for(int64_t p = 0 ; p < N ; p++)
    {
        // For each p, we substract Bg/2 * Bg^(-i) to pol_univ[p]
        tmp_pol_univ[p] = pol_univ[p];
        for(int64_t i = 1 ; i < l ; i++){
            tmp_pol_univ[p] += ldexp(1.0, kappa - 1 - kappa*i);
        }

        if(tmp_pol_univ[p] >= 0)
        {
            for(int64_t i = 0 ; i < l ; i++) //TODO to discuss
            {
                pol_biv[i*N + p] = (((int64_t)ldexp(tmp_pol_univ[p], i*kappa)) & mask) - (1LL << (kappa - 1));
            }
        }
        else{
            tmp_pol_univ[p] -= floor(tmp_pol_univ[p]);
            for(int64_t i = 0 ; i < l ; i++) //TODO to discuss
            {
                pol_biv[i*N + p] = (((int64_t)ldexp(tmp_pol_univ[p], i*kappa)) & mask) - (1LL << (kappa - 1));
            }
        }
    }

    free(tmp_pol_univ);

    return 0;
}