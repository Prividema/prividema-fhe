#include <string.h>
#include "math.h"
#include "ggsw.h"
#include "distributions.h" // Allow to draw random number following the uniform or normal law
#include "vec_znx_arithmetic_private.h"


/**
 * @brief Computes a random normal bivariate polynomial
 * 
 * @param module The module stocking the degree N.
 * @param kappa 
 * @param l 
 * @param res 
 * @param n_limbs 
 * @param res_sl 
 * 
 * @return int 
 */
int normal_random_biv_poly(MODULE* module, int64_t kappa, int64_t l, PolyBiv* res, int64_t n_limbs, int64_t res_sl
){
    int64_t N = module->nn;
    int64_t* tmp_biv_pol = malloc(N * l * sizeof(double));
    
    // Fills res(Y^0) with big coefficients
    for(int64_t p = 0 ; p < N ; p++)
    {
        if(rand_normal(((double *)tmp_biv_pol) + p, 0.0, 1.0) < 0) // TODO : decompose to xy
                return -1;

        tmp_biv_pol[p] = (int64_t) ldexp(tmp_biv_pol[p] - 0.5, 64);
    }

    uint8_t* tmp_space = vec_znx_normalize_base2k_tmp_bytes(module);

    // Then base-2K normalization 
    vec_znx_normalize_base2k(module, kappa, res, l, N, tmp_biv_pol, l, N, tmp_space);

    return 0;
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
 * @param module The module stocking the degree N.
 * @param res_dft The result in DFT space.
 * @param res_size The result's size.
 * @param c_dft The left-hand side polynomial in DFT space .
 * @param c_size The left-hand size of c_dft.
 * @param d_dft The right-hand side polynomial in DFT space.
 * @param d_size The right-hand size of c_dft.
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
 * @brief Decrypts the phase (message + noise) and puts it in phase.
 * 
 * @param phase The phase. 
 * @param key The secret key.
 * @param c The ciphertext.
 */
void decrypt_biv_glwe(int64_t* phase, 
                      GGSWSecretKey* key,
                      int64_t* c
){

}

/**
 * @brief For each `i in [0,l]`, Computes `∑_j{0,k-1}[DFT(sk_j * limb_i(a_j))]` and Adds it to `DFT(limb_i(b))`.
 * 
 * @param params GLWE parameters.
 * @param res_ct_dft The result ciphertext in DFT space.
 * @param sk_dft The secret key in DFT space.
 * 
 * @note `sk_dft = (DFT(sk_0), ... ,DFT(sk_{k-1}))`.
 */
void secret_key_mult_add_inplace_dft(GLWECtParams* params,
                         VecBivDFT* res_ct_dft,
                         VecUnivDFT* sk_dft
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    #ifndef WITHY0
    for (int i = 0; i < l + 1 ; i++)
    {
        // The pointer to (DFT(limb_i(a_0)) , DFT(limb_i(a_1)) , ... , DFT(limb_i(b))) 
        double* ct_limb_i_dft = res_ct_dft + i*N*(k+1);

        // Compute ∑_j{0,k-1}[DFT(sk_j * limb_i(a_j))] and Adds it to DFT(limb_i(a_k)) = DFT(limb_i(b))  
        for (int j = 0 ; j < k; j++)
        { 
            for (int64_t p = 0 ; p < N/2 ; j++)
            {
                
                // j*N + p corresponds to the p-th coefficient's Re[DFT(a_j_i)] & Re[DFT(sk_j)]
                double c_re = ct_limb_i_dft[j*N + p];
                double d_re = sk_dft[j*N + p];
                
                // j*N + p + N/2 corresponds to the p-th coefficient's Im[DFT(a_j_i)] & Im[DFT(sk_j)]
                double c_im = ct_limb_i_dft[j*N + p + N/2];
                double d_im = sk_dft[j*N + p + N/2];

                ct_limb_i_dft[k*N + j] += c_re*d_re - c_im*d_im;
                ct_limb_i_dft[k*N + j + N/2] += c_re*d_im + c_im*d_re;
            }
        }
    }
    #endif
}

/**
 * @brief Encrypts the phase (message + noise) and puts it in res.
 * 
 * @param module The module stocking the degree N.
 * @param params The GLWE parameters.
 * @param res_ct The result bivariate ciphertext. 
 * @param sk The secret key in DFT space.
 * @param phase message + noise.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */ // TODO false
int encrypt_biv_glwe(const MODULE* module,
                     GLWECtParams* params, 
                     VecBiv* res_ct,
                     GGSWPreparedSK* sk, 
                     PolyBiv* phase
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t kappa = params->kappa;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    if (uniform_random_vec(k * N, res_ct, l, (k + 1) * N) > 0 ) {
        return -1;
    }
    
    // acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    if (!acc){
        perror("calloc failed");
        return -1;
    }

    // Computes ∑_j{0,k-1}[resVec_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        SVP_PPOL* sk_j_dft = (SVP_PPOL *)sk->values[j]; 
        
        // Computes resVec_j_dft = (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j_dft = (PolyBivDFT *)new_vec_znx_dft(module, l); 
        svp_apply_dft(module, (VEC_ZNX_DFT *)resVec_j_dft, l, sk_j_dft, res_ct + j*N, l, (k+1)*N); 
        
        // Computes resVec_j in ZnXY space
        PolyBiv* resVec_j = (PolyBiv*)new_vec_znx_big(module, l); 
        uint8_t* tmp = 0x0;
        vec_znx_idft(module, (VEC_ZNX_BIG *)resVec_j, l, (VEC_ZNX_DFT *)resVec_j_dft, l, tmp);

        // And adds it to acc_j
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] += resVec_j[p];
        }
        delete_vec_znx_dft(resVec_j_dft);
        delete_vec_znx_big(resVec_j);
    }

    // The pointer to limb_0(b)
    PolyUnivDFT* b_0 = res_ct + k*N;
    uint8_t* tmp = vec_znx_normalize_base2k_tmp_bytes(module);

    // For each i in {0,l} limb_i(b) = acc_i = ∑_j{0,k-1}[s_j * limb_i(a_j)]
    vec_znx_normalize_base2k(module, kappa, b_0, l, N*(k+1), acc, l, N, tmp);
    
    free(tmp);
    free(acc);

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
                        GGSWPreparedSK* sk,             
                        PolyUniv* msg,                
                        GGSWCtParams* enc_params 
){
    // GLWE parameters
    int64_t N = res->params->params->N;
    int64_t n_limbs = res->params->params->n_limbs;

    // Matrix parameters
    //! The partials' structure is described in ggsw_ciphertext.dox
    int64_t nb_partial = res->params->n_limbs_tilde/(res->params->k_tilde + 1);
    int64_t nb_rows_per_partial = res->params->k_tilde + 1;
    
    MODULE* module = new_module_info(N,FFT64);
    
    // Computes message in DFT space
    PolyUnivDFT* msg_dft = (PolyUnivDFT*)new_vec_znx_dft(module, 1);
    vec_znx_dft(module, (VEC_ZNX_DFT*)msg_dft, 1, msg, 1, N);
    
    for (int64_t i = 0 ; i < nb_partial ; i++){
        for (int64_t j = 0 ; j < nb_rows_per_partial ; j++){

            // The pointer to bivGLWE(-m * s_j * Y^i)
            VecBiv* ct_biv = ggsw_Sj_Yi(res, i, j);
            
            // The pointer to DFT(sk_j)
            PolyUnivDFT* sk_j_dft = sk->values[j];
            
            // Compute -DFT(msg * sk_j)
            PolyUnivDFT* phase_dft = calloc(N,sizeof(double));
            if (phase_dft == NULL)
                return -1;
                
            vec_znx_dft_mult(module, (VEC_ZNX_DFT*)phase_dft, 1, (VEC_ZNX_DFT*)sk_j_dft, 1, (VEC_ZNX_DFT*)msg_dft, 1);
            for(int64_t p = 0 ; p < N ; p++){
                phase_dft[p] = -1 * phase_dft[p];  
            }

            // Compute the base-2Kappa of the phase -DFT(msg * sk_j)*Y^i
            PolyBiv* phase_biv = malloc(N * l * sizeof(int64_t));
            aux_compute_phase_biv(phase_biv); // TODO define it

            #ifdef WITH_Y0 
            if (encrypt_biv_glwe(module, res->params->params, 
                                 ct_biv, sk, phase_biv) < 0){
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

    delete_module_info(module);
    delete_vec_znx_dft(msg_dft);   
}

/**
 * @brief Encrypts the phase (message + noise) in DFT space and puts it in res_ct.
 * 
 * @param params The GLWE parameters.
 * @param res_ct The result ciphertext in DFT space. 
 * @param module The module stocking the degree N.
 * @param sk The secret key in DFT space.
 * @param phase message + noise in DFT space.
 * 
 * @retval - `-1` if an error occurs. In this case the error is from a syscall and perror is called.
 * @retval - `0` otherwise.
 */
int encrypt_biv_glwe_dft(GLWECtParams* params, 
                         VecBivDFT* res_ct,
                         const MODULE* module, 
                         GGSWPreparedSK* sk, 
                         PolyBivDFT* phase_dft
){
    int64_t N = params->N;
    int64_t k = params->k;
    int64_t kappa = params->kappa;
    int64_t n_limbs = params->n_limbs;
    int64_t l = n_limbs / (k+1);

    // Temporary bivGLWE ciphertext 
    VecBiv* tmp_ct = malloc(N*l*(k+1)*sizeof(int64_t));

    if (uniform_random_vec(k * N, tmp_ct, l, (k + 1) * N) < 0 )
    {
        free(tmp_ct);
        return -1;
    }
    
    // TODO coeff between 0 and Bg
    // acc_(j+1) = acc_j + (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    
    // Computes ∑_j{0,k-1}[resVec_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        SVP_PPOL* sk_j_dft = (SVP_PPOL*)sk->values[j]; 
        
        // Computes resVec_j_dft = (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j_dft = (PolyBivDFT *)new_vec_znx_dft(module, l); 
        svp_apply_dft(module, (VEC_ZNX_DFT*)resVec_j_dft, l, sk_j_dft, tmp_ct + j*N, l, (k+1)*N); 
        
        // Computes resVec_j in ZnXY space
        PolyBiv* resVec_j = (PolyBiv*)new_vec_znx_big(module, l); 
        uint8_t* tmp = 0x0;
        vec_znx_idft(module, (VEC_ZNX_BIG *)resVec_j, l, (VEC_ZNX_DFT *)resVec_j_dft, l, tmp);

        // And adds it to acc_j : acc_(j+1) = acc_j + resVec_j
        for(int64_t p = 0 ; p < N*l ; p++)
        {
            acc[p] += resVec_j[p];
        }
        delete_vec_znx_dft(resVec_j_dft);
        delete_vec_znx_big(resVec_j);
    }
    
    // The pointer to limb_0(b) in ZnXY
    PolyUnivDFT* b_0 = tmp_ct + k*N;
    uint8_t* tmp = vec_znx_normalize_base2k_tmp_bytes(module);

    // For each i in {0,l} limb_i(b) = acc_i = ∑_j{0,k-1}[s_j * limb_i(a_j)]
    // Then b is normalized
    vec_znx_normalize_base2k(module, kappa, b_0, l, N*(k+1), acc, l, N, tmp);
    
    // Computes tmp_ct in DFT space
    vec_znx_dft(module, (VEC_ZNX_DFT *)res_ct, l*(k+1), tmp_ct, l*(k+1), N);

    // Adds error to the phase
    add_error();
    for(int i = 0 ; i < l ; i++){
        PolyUniv* tmp = malloc(N*sizeof(double));
        PolyBiv* tmp_biv = malloc(N*l*sizeof(int64_t));

        if(normal_random_vec(N, tmp, 1, 1) < 0)
            return -1;

    }
    
    // Adds the phase (message with error) to bivGLWE(0), the result is a ct of bivGLWE(m + e)
    for (int64_t i = 0 ; i < l ; i++)
    {
        for (int64_t p = 0 ; p < N ; p++)
        {
            // Adds DFT(limb_i(phase)) to DFT(limb_i(b))
            res_ct[i*N*(k+1) + k*N + p] = res_ct[i*N*(k+1) + k*N + p] + phase_dft[i*N + p];
        }
    }
    
    free(acc);
    free(tmp_ct); 
    free(tmp);

    return 0;
}

/**
 * @brief Compute the base-2^kappa decomposition of the univariate polynomial in DFT space and return in DFT space
 * 
 * @param params The GLWE parameters.
 * @param pol_biv_dft The bivariate phase in DFT space
 * @param pol_dft The univariate phase in DFT space
 * @param i The degree in Y of DFT(-m * sk_j) * Y^i
 
void decomp_base_2K(GLWECtParams* params, PolyBivDFT* pol_biv_dft, PolyUnivDFT* pol_dft, int64_t precision){
*/


/**
 * @brief Compute the base-2^kappa decomposition of the univariate polynomial in DFT space and return in DFT space
 * 
 * @param params The GLWE parameters.
 * @param pol_biv_dft The bivariate phase in DFT space
 * @param pol_dft The univariate phase in DFT space
 * @param i The degree in Y of DFT(-m * sk_j) * Y^i
 
void decomp_base_2K_dft_to_dft(GLWECtParams* params, PolyBivDFT* pol_biv_dft, PolyUnivDFT* pol_dft, int64_t i){

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
*/

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
    //! The partials' structure is described in ggsw_ciphertext.dox
    int64_t nb_partial = n_limbs_tilde/(k_tilde + 1);
    int64_t nb_rows_per_partial = k_tilde + 1;
    
    // Prepare sk and m
    MODULE* module = new_module_info(N,FFT64);
    
    // Message, univariate polynomial in ZnX
    VEC_ZNX_DFT* msg_dft = new_vec_znx_dft(module, 1);
    vec_znx_dft(module, msg_dft, 1, msg, 1, N);

    for (int64_t i = 0 ; i < nb_partial ; i++)
    {
        for (int64_t j = 0 ; j < nb_rows_per_partial ; j++)
        {
            // The pointer to bivGLWE(-m * s_j * Y^i) in DFT space
            VecBivDFT* ct_biv_dft = res->ct + i*nb_rows_per_partial*n_limbs*N + j*n_limbs*N;
            
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

            delete_vec_znx_dft((VEC_ZNX_DFT*)phase_dft);
            
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