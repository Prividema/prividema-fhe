#include <string.h>
#include "math.h"
#include "ggsw.h"
#include "distributions.h" // Allow to draw random number following the uniform or normal law
#include "glwe_ciphertext.h"
#include "bivariate_polynomial.h"
#include "vec_znx_arithmetic_private.h"

//! GGSW PART (begin)

/**
 * @brief Decrypts the phase (message + noise) and puts it in phase.
 * 
 * @param enc_params The GLWE parameters.
 * @param phase The phase in RnX. 
 * @param key The secret key in DFT space.
 * @param ct The ciphertext.
 */
int decrypt_biv_glwe(GLWECtParams* enc_params,
                     double* phase, 
                     GGSWPreparedSK* sk_dft,
                     GLWECiphertext* ct
){
    // GLWE parameters
    int64_t N = enc_params->N;
    int64_t k = enc_params->k;
    int64_t l = poly_biv_size(enc_params);

    MODULE* module = new_module_info(N, FFT64);

    PolyBiv* acc = malloc(poly_biv_bytes(enc_params)); 
    if (!acc){
        perror("calloc failed");
        return -1;
    }

    // Computes acc = -∑_j{0,k-1}[sk_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of resp. the secret key and the bivGLWE ciphertext 
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        PolyUniv* a_j = ct->vec + j*N;
        
        // Computes DFT(sk_j * a_j)
        PolyBivDFT* as_j_dft = new_vec_znx_dft_p(module, l); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k+1)*N); 
        
        // Computes sk_j * a_j
        PolyBiv* as_j = new_vec_znx_big_p(module, l); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // And subs it to acc
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] -= as_j[p];
        }
        delete_vec_znx_dft_p(as_j_dft);
        delete_vec_znx_big_p(as_j);
    }

    // Computes acc = b - ∑_j{0,k-1}[sk_j * a_j]
    int64_t* b = ct->vec + N*k;
    add_biv_poly(enc_params, acc, N, b, N*(k+1), acc, N);
    
    biv_to_univ(enc_params, phase, acc);

    return 0;
}

void ggsw_decrypt(double* res,   // result
                  GGSWPreparedSK* sk_dft,  // secret key
                  GGSWCiphertext* ct  // ciphertext
){

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
                     GGSWPreparedSK* sk_dft, 
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

    // Computes ∑_j{0,k-1}[resVec_j_dft]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes resVec_j_dft = (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j_dft = new_vec_znx_dft_p(module, l); 
        svp_apply_dft_p(module, resVec_j_dft, l, sk_j_univ_dft, res_ct + j*N, l, (k+1)*N); 
        
        // Computes resVec_j in ZnXY space
        PolyBiv* resVec_j = new_vec_znx_big_p(module, l); 
        vec_znx_idft_p(module, resVec_j, l, resVec_j_dft, l);

        // And adds it to acc_j
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] += resVec_j[p];
        }
        delete_vec_znx_dft_p(resVec_j_dft);
        delete_vec_znx_big_p(resVec_j);
    }

    // The pointer to limb_0(b)
    PolyUniv* b_0_univ = res_ct + k*N;

    // For each i in {0,l} limb_i(b) = acc_i = ∑_j{0,k-1}[s_j * limb_i(a_j)]
    vec_znx_normalize_base2k_p(module, kappa, b_0_univ, l, N*(k+1), acc, l, N);
    
    free(acc);

    return 0;
}

/**
 * @brief Compute the base-2Kappa decomposition of the phase = -m * sk_j * (1/Bg_tilde)^i and put it in DFT space.
 * 
 * @param params The GGSW parameters.
 * @param res_biv_dft The result bivariate phase in DFT space.
 * @param phase_dft The initial univariate phase in DFT space.
 * @param i The exponent of 1/Bg_tilde.
 * @return int 
 * 
 * @note We look for the minimal p such that, p*kappa > i*kappa_tilde,
 * ie p = qot(i*kappa_tilde,kappa) + 1 and we compute `(-m * sk_j) * (Bg^p)/(Bg_tilde^i)`.
 * Finally we put it in front of Y^p and normalize.
 * 
 * If p > l_tilde then the phase is to low and we round it to zero.
 */
int compute_phase_biv(GGSWCtParams* params, 
                      PolyBiv* res, 
                      PolyUnivDFT* phase_univ_dft, 
                      int64_t i
){
    // GGSW parameters
    int64_t kappa_tilde = params->kappa_tilde;
    int64_t l_tilde = params->n_limbs_tilde/(params->k_tilde + 1);

    // GLWE parameters
    int64_t N = params->params->N;
    int64_t kappa = params->params->kappa;
    int64_t l = poly_biv_size(params->params);

    PolyUniv* phase_univ = malloc(poly_univ_bytes(params->params));
    if (phase_univ == NULL)
        perror("malloc failed");
        return -1;
    PolyBiv* phase_biv = calloc(poly_biv_bytes(params->params),1);
    if (phase_biv == NULL)
        perror("malloc_faioled");
        return -1;

    int64_t p = (i * kappa)/kappa_tilde + 1;

    if (p <= l_tilde)
    {
        MODULE* module = new_module_info(N, FFT64);

        // Compute iDFT(DFT(-m * sk_j)) = -m * sk_j
        vec_znx_idft_p(module, phase_univ, 1, phase_univ_dft, 1);
        
        // Compute -m * sk_j * (Bg_tilde^p)/(Bg^i) = -m * sk_j * 2^(p*kappa_tilde - i*kappa).
        // Then, put it in phase_biv at Y^p.
        for(int64_t t = 0 ; t < N ; t++) // TODO to discuss 
        {
            phase_univ[t] = phase_univ[t] << (p*kappa_tilde - i*kappa); 
            phase_biv[p*N + t] = phase_univ[t];
        }
        
        // Normalize phase_biv.
        
        vec_znx_normalize_base2k_p(module, kappa_tilde, res, l, N, phase_biv, l, N);

        free(module);
        free(phase_univ);
        free(phase_biv);
    }

    free(phase_univ);
    free(phase_biv);

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
                        GGSWPreparedSK* sk_dft,             
                        PolyUniv* msg_univ,                
                        GGSWCtParams* enc_params 
){
    // GGSW & GLWE parameters
    GGSWCtParams* params_ggsw = enc_params;
    GLWECtParams* params_glwe = enc_params->params;
    int64_t N = params_glwe->N;
    
    MODULE* module = new_module_info(N,FFT64);
    
    // Computes message in DFT space
    PolyUnivDFT* msg_univ_dft = new_vec_znx_dft_p(module, 1);
    vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, N);
    
    for (int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++){
        for (int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++){

            // The pointer to bivGLWE(-m * s_j * Y^i)
            VecBiv* ct_biv = ggsw_Sj_Yi(res, i, j);
            
            // The pointer to DFT(sk_j)
            PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
            
            // Compute -DFT(msg * sk_j)
            PolyUnivDFT* phase_univ_dft = calloc(N,sizeof(double));
            if (phase_univ_dft == NULL)
                return -1;
                
            vec_znx_dft_mult(module, phase_univ_dft, 1, sk_j_univ_dft, 1, msg_univ_dft, 1);
            for(int64_t p = 0 ; p < N ; p++){
                phase_univ_dft[p] = -1 * phase_univ_dft[p];  
            }
            
            // Compute the base-2^kappa decomposition of the phase = DFT(-m * sk_j) * (1/Bg)^i 
            // and return in DFT space
            PolyBiv* phase_biv = malloc(poly_biv_size(params_glwe) * sizeof(int64_t));
            compute_phase_biv(res->params, phase_biv, phase_univ_dft, i);

            #ifdef WITH_Y0 
            if (encrypt_biv_glwe(module, res->params->params, 
                                 ct_biv, sk_dft, phase_biv) < 0){
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
    delete_vec_znx_dft_p(msg_univ_dft);   

    return 0;
}


//! GGSW IN DFT PART (begin)

/**
 * @brief Adds a bivariate error to the bivariate phase, returns in DFT space.
 * 
 * @param enc_params The GLWE parameters.
 * @param res_dft The result bivariate phase in DFT space.
 * @param phase_dft The input phase.
 * @return int 
 */
int add_error_dft(GLWECtParams* enc_params,
                  PolyBivDFT* res_dft,
                  PolyBivDFT* phase_dft 
){
    MODULE* module = new_module_info(enc_params->N, FFT64);
    
    // Compute a random error in DFT space
    PolyBivDFT* err_dft = new_normal_random_biv_poly_dft(module, enc_params);
    
    // Add the error in DFT space
    add_biv_poly_dft(enc_params, phase_dft, enc_params->N, phase_dft, enc_params->N, err_dft, enc_params->N);
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
int encrypt_biv_glwe_dft(GLWECtParams* enc_params, 
                         const MODULE* module, 
                         VecBivDFT* res_dft,
                         GGSWPreparedSK* sk_dft, 
                         PolyBivDFT* phase_dft
){
    int64_t N = enc_params->N;
    int64_t k = enc_params->k;
    int64_t kappa = enc_params->kappa;
    int64_t n_limbs = enc_params->n_limbs;
    int64_t l = n_limbs / (k+1);

    // Temporary bivGLWE ciphertext 
    VecBiv* tmp_ct = malloc(N*l*(k+1)*sizeof(int64_t));
    
    // TODO coeff between 0 and Bg
    if (uniform_random_vec(k * N, tmp_ct, l, (k + 1) * N) < 0 )
    {
        free(tmp_ct);
        return -1;
    }
    
    // acc_(j+1) = acc_j + (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    
    // Computes ∑_j{0,k-1}[resVec_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes resVec_j_dft = (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j_dft = new_vec_znx_dft_p(module, l); 
        svp_apply_dft_p(module, resVec_j_dft, l, sk_j_univ_dft, tmp_ct + j*N, l, (k+1)*N); 
        
        // Computes resVec_j in ZnXY space
        PolyBiv* resVec_j = new_vec_znx_big_p(module, l); 
        vec_znx_idft_p(module, resVec_j, l, resVec_j_dft, l);

        // And adds it to acc_j : acc_(j+1) = acc_j + resVec_j
        for(int64_t p = 0 ; p < N*l ; p++)
        {
            acc[p] += resVec_j[p];
        }
        delete_vec_znx_dft_p(resVec_j_dft);
        delete_vec_znx_big_p(resVec_j);
    }
    
    // The pointer to limb_0(b) in ZnXY
    PolyUniv* b_0_univ = tmp_ct + k*N;

    // For each i in {0,l} limb_i(b) = acc_i = ∑_j{0,k-1}[s_j * limb_i(a_j)]
    // Then b is normalized
    vec_znx_normalize_base2k_p(module, kappa, b_0_univ, l, N*(k+1), acc, l, N);
    
    // Computes tmp_ct in DFT space
    vec_znx_dft_p(module, res_dft, l*(k+1), tmp_ct, l*(k+1), N);

    // Adds error to the phase
    add_error_dft(enc_params, phase_dft, phase_dft);
    
    // Adds the phase (message with error) to bivGLWE(0), the result is a ct of bivGLWE(m + e)
    for (int64_t i = 0 ; i < l ; i++)
    {
        for (int64_t p = 0 ; p < N ; p++)
        {
            // Adds DFT(limb_i(phase)) to DFT(limb_i(b))
            res_dft[i*N*(k+1) + k*N + p] = res_dft[i*N*(k+1) + k*N + p] + phase_dft[i*N + p];
        }
    }
    
    free(acc);
    free(tmp_ct); 

    return 0;
}


/**
 * @brief Compute the base-2Kappa decomposition of the phase = -m * sk_j * (1/Bg_tilde)^i and put it in DFT space.
 * 
 * @param params The GGSW parameters.
 * @param res_biv_dft The result bivariate phase in DFT space.
 * @param phase_dft The initial univariate phase in DFT space.
 * @param i The exponent of 1/Bg_tilde.
 * @return int 
 * 
 * @note We look for the minimal p such that, p*kappa > i*kappa_tilde,
 * ie p = qot(i*kappa_tilde,kappa) + 1 and we compute `(-m * sk_j) * (Bg^p)/(Bg_tilde^i)`.
 * Finally we put it in front of Y^p and normalize.
 * 
 * If p > l_tilde then the phase is to low and we round it to zero.
 */
int compute_phase_biv_dft(GGSWCtParams* enc_params, 
                          PolyBivDFT* res_dft, 
                          PolyUnivDFT* phase_univ_dft, 
                          int64_t i
){
    // GGSW parameters
    int64_t kappa_tilde = enc_params->kappa_tilde;
    int64_t l_tilde = enc_params->n_limbs_tilde/(enc_params->k_tilde + 1);

    // GLWE parameters
    int64_t N = enc_params->params->N;
    int64_t kappa = enc_params->params->kappa;
    int64_t l = poly_biv_size(enc_params->params);

    PolyUniv* phase_univ = malloc(poly_univ_bytes(enc_params->params));
    if (phase_univ == NULL)
        perror("malloc failed");
        return -1;
    PolyBiv* phase = calloc(poly_biv_bytes(enc_params->params),1);
    if (phase == NULL)
        perror("malloc_faioled");
        return -1;

    int64_t p = (i * kappa)/kappa_tilde + 1;

    if (p <= l_tilde)
    {
        MODULE* module = new_module_info(N, FFT64);

        // Compute iDFT(DFT(-m * sk_j)) = -m * sk_j
        vec_znx_idft_p(module, phase_univ, 1, phase_univ_dft, 1);
        
        // Compute -m * sk_j * (Bg_tilde^p)/(Bg^i) = -m * sk_j * 2^(p*kappa_tilde - i*kappa).
        // Then, put it in phase_biv at Y^p.
        for(int64_t t = 0 ; t < N ; t++) // TODO to discuss 
        {
            phase_univ[t] = phase_univ[t] << (p*kappa_tilde - i*kappa); 
            phase[p*N + t] = phase_univ[t];
        }
        
        // Normalize phase_biv.
        
        vec_znx_normalize_base2k_p(module, kappa_tilde, phase, l, N, phase, l, N);
        vec_znx_dft_p(module, res_dft, l, phase, l, N);

        free(module);
        free(phase_univ);
        free(phase);
    }

    free(phase);
    free(phase_univ);

    return 0;

}

int ggsw_secret_encrypt_dft(GGSWCtParams* enc_params,    // parameters
                            GGSWCiphertextDFT* res_dft,        // result
                            GGSWPreparedSK* sk_dft,         // secret key
                            PolyUniv* msg_univ                 // message
){
    //def a
    //(a, sk *a ) + mu*ID n_limbs_tilde fois 

    // GGSW & GLWE parameters
    GGSWCtParams* params_ggsw = enc_params;
    GLWECtParams* params_glwe = enc_params->params;
    int64_t N = params_glwe->N;

    // Prepare sk and m
    MODULE* module = new_module_info(N,FFT64);
    
    // Message, univariate polynomial in ZnX
    PolyUnivDFT* msg_univ_dft = new_vec_znx_dft_p(module, 1);
    vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, N);

    for (int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++)
    {
        for (int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++)
        {
            // The pointer to bivGLWE(-m * s_j * Y^i) in DFT space
            VecBivDFT* ct_biv_dft = ggsw_Sj_Yi_dft(res_dft, i, j);
            
            // The pointer to DFT(sk_j)
            PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
            
            // Computes DFT(-m * sk_j)
            PolyUnivDFT* phase_univ_dft = calloc(N,sizeof(double));
            vec_znx_dft_mult(module, phase_univ_dft, 1, sk_j_univ_dft, 1, msg_univ_dft, 1);

            PolyBivDFT* phase_dft = malloc(poly_biv_bytes(params_glwe));

            // Compute the base-2^kappa decomposition of the phase = DFT(-m * sk_j) * (1/Bg)^i 
            // and return in DFT space
            compute_phase_biv_dft(res_dft->params, phase_dft, phase_dft, i);
            
            #ifdef WITH_Y0 
            // Compute bivGLWE(-m * s_j * Y^i)
            if (encrypt_biv_glwe_dft(res_dft->params->params, module, ct_biv_dft,
                                  sk_dft, phase_dft) < 0)
            {
                delete_vec_znx_dft_p(msg_univ_dft);
                delete_module_info(module);
                return -1;
            }

            delete_vec_znx_dft_p(phase_dft);
            free(phase_dft);
            
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
    delete_vec_znx_dft_p(msg_univ_dft);
    delete_module_info(module);

    return 0;
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