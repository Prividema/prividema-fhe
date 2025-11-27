#include "math.h"
#include "ggsw.h"
#include "distributions.h" // Allow to draw random number following the uniform or normal law

/**
 * @brief Compute the sum of two GGSW ciphertext
 * 
 * @param res the sum   
 * @param c first term of the sum
 * @param d second term of the sum
 */
void add_vec_znxy(GGSWCiphertext* res, GGSWCiphertext* c, GGSWCiphertext* d){

    int64_t nb_rows = res->params->n_limbs_tilde;
    int64_t nb_cols = res->params->params->n_limbs;
    int64_t N = res->params->params->N;

    for (int64_t i = 0 ; i < nb_rows ; i++){
        for (int64_t j = 0 ; j < nb_cols ; j++){
            for (int64_t k = 0 ; k < N ; k++){
                res->ct[i*N*nb_cols + j*N + k] = c->ct[i*N*nb_cols + j*N + k] + d->ct[i*N*nb_cols + j*N + k];
            } 
        }
    }
}

void add_vec_znx(int64_t* res, int64_t res_size, 
                 int64_t* c  , int64_t c_size  ,
                 int64_t* d  , int64_t d_size
){
    
    if (d_size > c_size){
        int64_t min_size = res_size > c_size ? c_size : res_size;
        int64_t add_size = res_size > d_size ? d_size : res_size;
        for (int64_t i = 0 ; i < add_size ; i++){
            res[i] = c[i] + d[i];
        }
        for (int64_t i = add_size ; i < min_size ; i++){
            res[i] = c[i];
        }
        for (int64_t i = min_size ; i < res_size ; i++){
            res[i] = 0;
        }
    }
    else {
        add_vec_znx(res, res_size, d, d_size, c, c_size);
    }
}

/**
 * @brief Encrypts the phase (message + noise) and puts it in res
 * 
 * @param res The result
 * @param sk The secret key
 * @param phase message + noise
 * @param encrypt_zero 1 if enrypting zero, 0 otherwise
 */
void encrypt_biv_glwe(int64_t* res, 
                      GGSWSecretKey* sk, 
                      int64_t* phase,
                      int encrypt_zero
){

}

/**
 * @brief Encrypts message m into GGSW ciphertext res with parameters enc_params
 * 
 * @param res The encrypted message
 * @param sk The secret key
 * @param m The message
 * @param enc_params The encryption params
 */
void ggsw_secret_encrypt(GGSWCiphertext* res,           
                         GGSWSecretKey* sk,             
                         int64_t* m,                
                         GGSWCtParams* enc_params 
){
    //def a
    //(a, sk *a ) + mu*ID n_limbs_tilde fois 

    // GLWE parameters
    int64_t N = res->params->params->N;
    int64_t k = res->params->params->k;
    int64_t n_limbs = res->params->params->n_limbs;

    // GGSW parameters
    int64_t k_tilde = res->params->k_tilde;
    int64_t n_limbs_tilde = res->params->n_limbs_tilde;

    // Matrix parameters
    int64_t nb_partial = k_tilde;
    int64_t nb_rows_per_partial = n_limbs_tilde/k_tilde;

    for (int64_t i = 0 ; i < nb_partial ; i++){
        for (int64_t j = 0 ; j < nb_rows_per_partial ; j++){
            int64_t* mm;
            #ifdef WITH_Y0 
            encrypt_biv_glwe(res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs, sk, mm, 0);
            add_vec_znx(res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs + i*(k+1)*N + j*N, N, 
                        res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs + i*(k+1)*N + j*N, N, 
                        m, N);
            #endif 
            #ifndef WITH_Y0
            encrypt_biv_glwe_without_y0(res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs, sk, mm, 0);
            add_vec_znx(res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs + (i - 1)*(k+1)*N + j*N, N, 
                        res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs + (i - 1)*(k+1)*N + j*N, N, 
                        m, N);
            #endif 
        }
    }

}

int* add(int* a, int a_size, int* b, int b_size) {
    if (a_size >= b_size){
        for (int i = 0; i < b_size ; i++){
            a[i] = a[i] + b[i];
        }
        return a;
    }
    else {
        return add(b,b_size,a,a_size);
    }
}

int add_random_int(int a, int b){
    return a+b;}
int multiply(int a, int b) {
    return a * b;
}