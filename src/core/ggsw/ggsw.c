#include <string.h>
#include "math.h"
#include "ggsw.h"
#include "distributions.h" // Allow to draw random number following the uniform or normal law
#include "vec_znx_arithmetic_private.h"

/* Adds two GGSW ciphertexts with same params and put result in res */
void ggsw_add(GGSWCiphertext* res,  // result
             GGSWCiphertext* ct1,  // first operand
             GGSWCiphertext* ct2   // second operand
){
    int64_t nb_rows = res->params->n_limbs_tilde;
    int64_t nb_cols = res->params->params->n_limbs;
    int64_t N = res->params->params->N;

    for (int64_t i = 0 ; i < nb_rows ; i++){
        for (int64_t j = 0 ; j < nb_cols ; j++){
            for (int64_t k = 0 ; k < N ; k++){
                res->ct[i*N*nb_cols + j*N + k] = ct1->ct[i*N*nb_cols + j*N + k] + ct2->ct[i*N*nb_cols + j*N + k];
            } 
        }
    }
}


/**
 * @brief Adds the bivariate polynomial m to the GGSW bivariate ciphertext 
 * 
 * @details If ct_in = (a_0,..,a_p,..,a_k-1,b) then ct_in = (a_0,..,a_p + m,..,a_k-1,b)
 * 
 * @param params GGSW parameters
 * @param ct_out GGSW bivariate ciphertext output
 * @param m bivariate polynomial
 * @param ct_in GGSW bivariate ciphertext input
 * @param p the index of `a_p`
 * @note If p = k then `a_k = b` 
 *
 */

void add_inplace_m_to_aj_dft(GLWECtParams* params, 
                         VecBivDFT* ct_in,
                         PolyBivDFT* m,
                         int64_t j
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t l = params->n_limbs/(k+1); 

    for (int64_t i = 0 ; i < l ; i++){
        for (int64_t p = 0 ; p < N ; p++){
            // No need to differenciate Re and Im since we are only adding m to ct_in
            ct_in[i*N*(k+1) + j*N + p] = ct_in[i*N*(k+1) + j*N + p] + m[j*N + p];
        }
    }
}

// TODO define multiplication for each type
/**
 * @brief Compute the polynomial product of c and d, component-wise in DFT space.
 * 
 * @param module The module holding N the degree in X and 
 * @param res_dft The result in DFT space
 * @param res_size The result's size
 * @param c_dft The left-hand side polynomial in DFT space 
 * @param c_size 
 * @param d_dft 
 * @param d_size 
 * 
 * @note `res_dft = ( DFT(c_0) * DFT(d_0) , ... , DFT(c_smin) * DFT(d_smin) , 0's)`. There are enough 0's to match the size of res_dft.
 */
void vec_znx_dft_mult(const MODULE* module, 
              VEC_ZNX_DFT* res_dft, int64_t res_size,
              VEC_ZNX_DFT* c_dft, int64_t c_size,  
              VEC_ZNX_DFT* d_dft, int64_t d_size
){
    int64_t N = module->nn;

    if (c_size < d_size){
        int64_t smin = c_size < res_size ? c_size : res_size;
        
        for (int i = 0 ; i < smin; i++){
            for (int64_t j = 0 ; j < N/2 ; j++){ // TODO
                // i*N + j corresponds to the p-th coefficient's Im[DFT(c)] & Im[DFT(d)]
                double c_re = ((double*)c_dft)[i*N + j];
                double d_re = ((double*)d_dft)[i*N + j];

                // i*N + j + N/2 corresponds to the p-th coefficient's Im[DFT(c)] & Im[DFT(d)]
                double c_im = ((double *)c_dft)[i*N + j + N/2];
                double d_im = ((double *)d_dft)[i*N + j];

                ((double*)res_dft)[i*N + j] = c_re * d_re - c_im * d_im;
                ((double*)res_dft)[i*N + j + N/2] = c_re * d_im + c_im * d_re; 
            }
        }
        
        // fill up remaining part with 0's
        double* const dres_dft = (double*)res_dft;
        memset(dres_dft + smin* N, 0, (res_size - smin) * N * sizeof(double));
    }
    else {
        vec_znx_dft_mult(module, res_dft, res_size, d_dft, d_size, c_dft, c_size);
    }
}

/**
 * @brief Compute `∑_j{0,k-1}[ sk_j * a_j_i ]` and Adds it to `a_k_i = b_i` in place, in DFT space
 * 
 * @param module The module holding N the degree in X and 
 * @param a_j0k_i_dft The GLWE ciphertext's `Y^i` precision component in DFT space
 * @param sk_dft The secret key in DFT space
 * 
 * @note `sk_dft = ( DFT(sk_0) , ... , DFT(sk_k) )`
 * @note `a_j0k_i_dft = ( DFT(a_0_i) , ... , DFT(a_k_i))`
 * 
 * @details ` DFT(b_i) = ∑_j{0,k-1}[ DFT(sk_j) * DFT(a_j_i) ]`
 */
void dft_mult_add_j0k_i_inplace(GLWECtParams* params,
                        VecUnivDFT* a_j0k_i_dft,  
                        VecUnivDFT* sk_dft
){  
    int64_t N = params->N;
    int64_t k = params->k;

    for (int j = 0 ; j < k; j++){ 
        for (int64_t p = 0 ; p < N/2 ; j++){
            
            // j*N + p corresponds to the p-th coefficient's Re[DFT(a_j_i)] & Re[DFT(sk_j)]
            double c_re = a_j0k_i_dft[j*N + p];
            double d_re = sk_dft[j*N + p];
            
            // j*N + p + N/2 corresponds to the p-th coefficient's Im[DFT(a_j_i)] & Im[DFT(sk_j)]
            double c_im = a_j0k_i_dft[j*N + p + N/2];
            double d_im = sk_dft[j*N + p + N/2];

            a_j0k_i_dft[k*N + j] += c_re*d_re - c_im*d_im;
            a_j0k_i_dft[k*N + j + N/2] += c_re*d_im + c_im*d_re;
        
        }
    }
}

/**
 * @brief Decrypts the phase (message + noise) and puts it in phase
 * 
 * @param phase The phase 
 * @param key The secret key
 * @param c The ciphertext
 */
void decrypt_biv_glwe(int64_t* phase, 
                      GGSWSecretKey* key,
                      int64_t* c
){

}

/**
 * @brief For each `i in [0,l]`, Computes `∑_j{0,k-1}[ DFT(sk_j * a_j_i) ]` and Adds it to `DFT(a_k_i) = DFT(b_i)`
 * 
 * @param params GLWE params 
 * @param res_ct_dft The result ciphertext in DFT space
 * @param sk_dft The secret key in DFT space
 * 
 * @note `sk_dft = (DFT(sk_0), ... ,DFT(sk_k))`
 */
void secret_key_mult_add_inplace(GLWECtParams* params,
                         double* res_ct_dft,
                         double* sk_dft
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    #ifndef WITHY0
    for (int j = 0; j < l + 1 ; j++){

        // The pointer to (DFT(a_0_i), DFT(a_1_i) , ... , DFT(b_i)) 
        double* a_j0k_i_dft = res_ct_dft + j*N*(k+1);

        // Compute ∑_j{0,k-1}( DFT(sk_j * a_j_i) ) and Adds it to DFT(a_k_i) = DFT(b_i)  
        dft_mult_add_j0k_i_inplace(params, a_j0k_i_dft, sk_dft);
    }
    #endif
}

/**
 * @brief Encrypts the phase (message + noise) and puts it in res.
 * 
 * @param ct The result ciphertext. 
 * @param sk The secret key.
 * @param phase message + noise.
 * @param encrypt_zero 1 if enrypting zero, 0 otherwise.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int encrypt_biv_glwe(GLWECtParams* params, 
                     int64_t* res_ct,
                     const MODULE* module, 
                     VEC_ZNX_DFT* sk_dft, int64_t sk_size, 
                     VEC_ZNX_DFT* phase
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    if (uniform_random_vec(k * N, res_ct , l, (k + 1) * N) > 0 ) {
        return -1;
    }


    VecBivDFT* res_ct_dft = new_vec_znx_dft(module, (k+1)*l);
    vec_znx_dft(module, res_ct_dft, (k+1)*l, res_ct, (k+1)*l, N);

    // Computes 
    secret_key_mult_add_inplace(params, (double *)res_ct_dft, (double *)sk_dft);
    
    add_inplace_m_to_ap_dft(params, (double *)res_ct_dft, (double *)phase, k);

    uint8_t* tmp;
    vec_znx_idft(module, (VEC_ZNX_BIG*)res_ct, (k+1)*l, res_ct_dft, (k+1)*l, tmp);

    delete_vec_znx_dft(res_ct_dft);
    return 0;
}

/**
 * @brief Encrypts message m into GGSW ciphertext res with parameters enc_params
 * 
 * @param res The encrypted message
 * @param sk The secret key
 * @param m The message
 * @param enc_params The encryption params
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int ggsw_secret_encrypt(GGSWCiphertext* res,           
                        GGSWSecretKey* sk,             
                        PolyUniv* m,                
                        GGSWCtParams* enc_params 
){
    //def a
    //(a, sk *a ) + mu*ID n_limbs_tilde fois 

    // GLWE parameters
    int64_t N = res->params->params->N;
    int64_t k = res->params->params->k;
    int64_t n_limbs = res->params->params->n_limbs;
    int64_t l = n_limbs / (k+1);

    // GGSW parameters
    int64_t k_tilde = res->params->k_tilde;
    int64_t n_limbs_tilde = res->params->n_limbs_tilde;

    // Matrix parameters
    int64_t nb_partial = k_tilde;
    int64_t nb_rows_per_partial = n_limbs_tilde/k_tilde;
    
    // Prepare sk and m
    MODULE* module = new_module_info(N,FFT64);
    
    // Secret key
    VEC_ZNX_DFT* sk_dft = new_vec_znx_dft(module, n_limbs);
    vec_znx_dft(module, sk_dft, sk->size, sk->values, sk->size, N);

    // Message
    VEC_ZNX_DFT* m_dft = new_vec_znx_dft(module, 1);
    vec_znx_dft(module, m_dft, 1, m, 1, N);
    
    // m * sk_j
    VEC_ZNX_DFT* m_sk_j_dft = new_vec_znx_dft(module,1);
    
    for (int64_t i = 0 ; i < nb_partial ; i++){
        for (int64_t j = 0 ; j < nb_rows_per_partial ; j++){

            // The pointer to bivGLWE(-m * s_j * Y^i)
            int64_t* ct_m_sj_Yi = res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs;
            
            // The pointer to DFT(sk_j)
            VEC_ZNX_DFT* sk_j_dft = (VEC_ZNX_DFT*)((double*)(sk_dft) + j*N);
            
            // Compute -m * s_j
            vec_znx_dft_mult(module, m_sk_j_dft, 1, sk_j_dft, 1, m_dft, 1);
            for(int64_t p = 0 ; p < N ; p++){
                ((double*)m_sk_j_dft)[p] = -1 * ((double*)m_sk_j_dft)[p];  
            }

            #ifdef WITH_Y0 
            if (encrypt_biv_glwe(res->params->params, ct_m_sj_Yi,
                                 module, sk_dft, sk->size, 
                                 m_sk_j_dft) < 0){
                return -1;
            }
            
            #endif 
            #ifndef WITH_Y0
            if (encrypt_biv_glwe(res->params->params, ct_m_sj_Yi,
                                 module, sk_dft, sk->size, 
                                 m_sk_j_dft) < 0){
                return -1;
            }
            #endif 
        }
    }

    delete_vec_znx_dft(m_dft);
    delete_vec_znx_dft(sk_dft);
    delete_vec_znx_dft(m_sk_j_dft);
    delete_module_info(module);
}

/**
 * @brief Encrypts the phase (message + noise) in DFT space and puts it in res.
 * 
 * @param ct The result ciphertext in DFT space. 
 * @param sk The secret key in DFT space.
 * @param phase message + noise in DFT space.
 * @param encrypt_zero 1 if enrypting zero, 0 otherwise.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int encrypt_biv_glwe_dft(GLWECtParams* params, 
                         VecBivDFT* res_ct,
                         const MODULE* module, 
                         GGSWPreparedSK* sk, 
                         PolyBivDFT* phase
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    int64_t* temp = calloc(N*l*(k+1),sizeof(int64_t));

    if (uniform_random_vec(k * N, temp , l, (k + 1) * N) > 0 )
        return -1;
    
    // acc_(j+1) = acc_j + (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
    double* acc = malloc(N*l*sizeof(double)); 
    if (!acc){
        perror("malloc failed");
        return -1;
    }
    
    // Computes ∑_j{0,k-1}[resVec_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        SVP_PPOL* sk_j_dft = (SVP_PPOL*)sk->values[j]; 
        
        // Computes resVec_j = (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j = new_vec_znx_dft(module, l); 
        svp_apply_dft(module, resVec_j, l, sk_j_dft, temp + j*N, l, (k+1)*N); 
        
        // And adds it to limb_i(b), for each i in {0,l}
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] += ((double*)resVec_j)[p];
        }
        delete_vec_znx_dft(resVec_j);
    }
    
    // Computes res in DFT space
    vec_znx_dft(module, res_ct, l*(k+1), temp, l*(k+1), N);

    // For each i in {0,l} limb_i(b) = ∑_j{0,k-1}[DFT(s_j) * limb_i(a_j)] = acc_i
    for(int64_t i = 0 ; i < l ; i++)
    {
        // The pointer to limb_i(b)
        double* b_i = ((double*)res_ct) + (k+1)*N*i + k*N;
        for(int64_t p = 0 ; p < N ; p++)
        {
            b_i[p] = acc[i*N + p];
        }
    }

    // Adds error to the phase
    // TODO

    // Adds the phase (message with error) to bivGLWE(0), the result is a ct of bivGLWE(m + e)
    add_inplace_m_to_ap_dft(params, res_ct, phase, k);
    
    free(acc);
    free(temp); 

    return 0;
}

/**
 * @brief Compute the base-2^kappa decomposition of the univariate phase and return in DFT space
 * 
 * @param params The GLWE parameters.
 * @param phase_biv_dft The bivariate phase in DFT space
 * @param phase_dft The univariate phase in DFT space
 * @param i The degree in Y of DFT(-m * sk_j) * Y^i
 */
void aux(GLWECtParams* params, PolyBivDFT* phase_biv_dft, PolyUnivDFT* phase_dft, int64_t i){

    int64_t N = params->N;
    int64_t k = params->k;
    int64_t kappa = params->kappa;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    MODULE* module = new_module_info(N,FFT64);
    VEC_ZNX_BIG* phase = new_vec_znx_big(module, 1);
    uint8_t* useless;
    vec_znx_idft(module, phase, 1, phase_dft, 1, useless);
    
    // The bivariate phase
    int64_t* phase_biv = calloc(N,sizeof(int64_t));

    // Computes the bivariate phase
    for(int64_t p = 0 ; p < N; p++)
    {
        phase_biv[i*N + p] = -1 * (int64_t)((double*)phase)[p]; // TODO does it works to do ((int64_t*)phase)[p]
    }
    // TODO Computes the bivariate phase
    vec_znx_dft(module, phase_biv_dft, l, phase_biv, l, N);
}

int ggsw_secret_encrypt_dft(GGSWPreparedCt* res,        // result
                            GGSWPreparedSK* sk,         // secret key
                            PolyUniv* msg,                 // message
                            GGSWCtParams* enc_params    // parameters
){
    //def a
    //(a, sk *a ) + mu*ID n_limbs_tilde fois 

    // GLWE parameters
    int64_t N = res->params->params->N;
    int64_t k = res->params->params->k;
    int64_t n_limbs = res->params->params->n_limbs;
    int64_t l = n_limbs / (k+1);

    // GGSW parameters
    int64_t k_tilde = res->params->k_tilde;
    int64_t n_limbs_tilde = res->params->n_limbs_tilde;

    // Matrix parameters
    int64_t nb_partial = k_tilde;
    int64_t nb_rows_per_partial = n_limbs_tilde/k_tilde;
    
    // Prepare sk and m
    MODULE* module = new_module_info(N,FFT64);
    
    // Message
    VEC_ZNX_DFT* msg_dft = new_vec_znx_dft(module, 1);
    vec_znx_dft(module, msg_dft, 1, msg, 1, N);

    for (int64_t i = 0 ; i < nb_partial ; i++)
    {
        for (int64_t j = 0 ; j < nb_rows_per_partial ; j++)
        {
            // The pointer to bivGLWE(-m * s_j * Y^i) in DFT space
            VecBivDFT* ct_biv_dft = res->ct + i*nb_rows_per_partial*n_limbs + j*n_limbs;
            
            // The pointer to DFT(sk_j)
            PolyUnivDFT* sk_j_dft = sk->values[j];
            
            // Computes DFT(-m * sk_j) * Y^i
            PolyUnivDFT* phase_dft = calloc(N,sizeof(double));
            vec_znx_dft_mult(module, (VEC_ZNX_DFT*)phase_dft, 1, (VEC_ZNX_DFT*)sk_j_dft, 1, msg_dft, 1);

            PolyBivDFT* phase_biv_dft = calloc(N*l,sizeof(double));


            // TODO
            // Compute the base-2^kappa decomposition of the phase and return in DFT space
            aux(res->params->params, phase_biv_dft, phase_dft, i);

            #ifdef WITH_Y0 
            // Compute bivGLWE(-m * s_j * Y^i)
            if (encrypt_biv_glwe_dft(res->params->params, ct_biv_dft,
                                 module, sk->values, phase_biv_dft) < 0)
            {
                delete_vec_znx_dft(msg_dft);
                delete_module_info(module);
                return -1;
            }

            delete_vec_znx_dft(phase_dft);
            
            #endif 
            #ifndef WITH_Y0
            if (encrypt_biv_glwe(res->params->params, ct_m_sj_Yi,
                                 module, sk_dft, sk->size, 
                                 m_sk_j_dft) < 0)
            {
                delete_vec_znx_dft(m_dft);
                delete_vec_znx_dft(m_sk_j_dft);
                delete_module_info(module);
                return -1;
            }
            #endif 
        }
    }
    delete_vec_znx_dft(msg_dft);
    delete_module_info(module);
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