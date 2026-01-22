#include "glwe_ciphertext.h"
#include "spqlios_alias.h"
#include "vec_znx_arithmetic_private.h"
#include <string.h>


//! GLWE PART (begin)

int64_t glwe_coef_number(GLWECtParams* params){
    return glwe_size(params) * params->N;
}

GLWECiphertext* new_glwe(GLWECtParams* params
){
    GLWECiphertext* ct = malloc(sizeof(GLWECiphertext));
    if(ct == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    ct->params = params;

    ct->vec = calloc(glwe_bytes(params), 1);
    if(ct->vec == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    return ct;
}

void delete_glwe(GLWECiphertext* ct){
    free(ct->vec);
    free(ct);
}


//! GLWE IN DFT PART (begin)

int64_t glwe_coef_number_dft(GLWECtParams* params){
    return glwe_size(params) * params->N / 2;
}

GLWEPreparedCt* new_glwe_dft(GLWECtParams* params){
    GLWEPreparedCt* ct = malloc(sizeof(GLWEPreparedCt));
    if(ct == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    ct->params = params;

    ct->pvec = calloc(glwe_bytes(params), 1);
    if(ct->pvec == NULL){
        perror("Malloc failed.");
        return NULL;
    }

    return ct;
}

void delete_glwe_dft(GLWEPreparedCt* ct){
    free(ct->pvec);
    free(ct);
}


//! COMMON PART (begin)

int64_t glwe_size(GLWECtParams* params){
    return params->n_limbs;
}


int64_t glwe_bytes(GLWECtParams* params){
    int64_t N = params->N;
    return glwe_size(params) * N * sizeof(int64_t); 
}


void mult_vec_znx_dft(const MODULE* module, 
              double* res_dft, int64_t res_size,
              double* c_dft, int64_t c_size,  
              double* d_dft, int64_t d_size
){
    int64_t N = module->nn;

    if (c_size <= d_size){
        int64_t smin = c_size < res_size ? c_size : res_size;
        
        for (int i = 0 ; i < smin; i++){
            for (int64_t j = 0 ; j < N/2 ; j++){ 
                // i*N + j corresponds to the j-th coefficient's index of Re[DFT(c_i)] and Re[DFT(d_i)]
                double c_re = c_dft[i*N + j];
                double d_re = d_dft[i*N + j];

                // i*N + j + N/2 corresponds to the j-th coefficient's index of Im[DFT(c_i)] & Im[DFT(d_i)]
                double c_im = c_dft[i*N + j + N/2];
                double d_im = d_dft[i*N + j + N/2];

                res_dft[i*N + j] = c_re * d_re - c_im * d_im;
                res_dft[i*N + j + N/2] = c_re * d_im + c_im * d_re; 
            }
        }
        
        // fill up remaining part with 0's
        double* const dres_dft = (double*)res_dft;
        memset(dres_dft + smin* N, 0, (res_size - smin) * N * sizeof(double));
    }
    else {
        mult_vec_znx_dft(module, res_dft, res_size, d_dft, d_size, c_dft, c_size);
    }
}