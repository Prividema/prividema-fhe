#include <string.h>
#include "math.h"
#include "ggsw.h"
#include "rng.h"
#include "glwe_ciphertext.h"
#include "bivariate_polynomial.h"
#include "vec_znx_arithmetic_private.h"

//! GGSW PART (begin)

void add_error(GLWECtParams* enc_params,
               PolyBiv* res,
               PolyBiv* phase)
{
    MODULE* module = new_module_info(enc_params->N, FFT64);
    
    // Compute a random error in DFT space
    PolyBiv* err = new_normal_random_biv_poly(module, enc_params);
    
    // Add the error in DFT space
    add_biv_poly(enc_params, res, enc_params->N, phase, enc_params->N, err, enc_params->N);

    free(err);
    delete_module_info_p(module);
}

int glwe_secret_demasking_ggsw_lib(GLWECtParams* enc_params,
                          double* phase,
                          GGSWSecretKeyDFT* sk_dft,
                          VecBiv* ct)
{
    // GLWE parameters
    uint64_t N = enc_params->N;
    uint64_t k = enc_params->k;
    uint64_t l = poly_biv_size(enc_params);

    MODULE* module = new_module_info(N, FFT64);

    PolyBiv* acc = malloc(poly_biv_bytes(enc_params)); 
    if (!acc){
        perror("calloc failed");
        return -1;
    }

    // Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of resp. the secret key and the bivGLWE ciphertext 
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        PolyUniv* a_j = ct + j*N;
        
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

    // Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
    int64_t* b = ct + N*k;
    add_biv_poly(enc_params, acc, N, b, N*(k+1), acc, N);
    
    PolyBiv* acc_normalized = malloc(poly_biv_bytes(enc_params));
    vec_znx_normalize_base2k_p(module, enc_params->kappa, acc_normalized, l, N, acc, l, N);

    biv_to_univ(enc_params, phase, acc_normalized);
    
    free(acc); free(acc_normalized);
    delete_module_info(module);

    return 0;
}

/**
 * @brief Computes Sum_j{0,k-1}[s_j * a_j].
 * 
 * @param module 
 * @param params 
 * @param res 
 * @param as 
 * @param sk_dft 
 * @return int 
 */
// TODO
int add_mult_ggsw(const MODULE* module, GLWECtParams* params,
             PolyBiv* res, VecBiv* as, GGSWSecretKeyDFT* sk_dft
){
    uint64_t N = params->N;
    uint64_t k = params->k;
    uint64_t l = poly_biv_size(params);

    // Computes Sum_j{0,k-1}[s_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes DFT(s_j) * DFT(a_j)
        PolyBivDFT* as_j_dft = calloc(poly_biv_coef_number_dft(params), 2*sizeof(double)); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, as + j*N, l, (k+1)*N); 
        
        // Computes s_j * a_j
        PolyBiv* as_j = calloc(poly_biv_coef_number(params), sizeof(int64_t)); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // And adds it to acc_j
        for(int64_t p = 0 ; p < N*l ; p++){
            res[p] += as_j[p];
        }
        free(as_j_dft); free(as_j);
    }
}

// TODO false
int glwe_secret_masking_ggsw_lib(const MODULE* module,
                        GLWECtParams* params, 
                        VecBiv* res_ct,
                        GGSWSecretKeyDFT* sk_dft, 
                        PolyBiv* phase)
{
    uint64_t N = params->N;
    uint64_t k = params->k;
    uint64_t kappa = params->kappa;
    uint64_t n_limbs = params->n_limbs;
    uint64_t l = n_limbs / (k+1);

    if (inplace_uniform_random_vec(k * N, res_ct, l, (k + 1) * N, kappa) > 0)
        return -1;
    
    // acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    if (!acc){
        perror("calloc failed");
        return -1;
    }

    // Computes Sum_j{0,k-1}[s_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes DFT(s_j) * DFT(a_j)
        PolyBivDFT* as_j_dft = calloc(poly_biv_coef_number_dft(params), 2*sizeof(double)); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, res_ct + j*N, l, (k+1)*N); 
        
        // Computes s_j * a_j
        PolyBiv* as_j = calloc(poly_biv_coef_number(params), sizeof(int64_t)); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // And adds it to acc_j
        for(int64_t p = 0 ; p < N*l ; p++){
            acc[p] += as_j[p];
        }
        free(as_j_dft); free(as_j);
    }

    // Add the phase to acc
    for(int64_t i = 0 ; i < l ; i++){    
        for(int64_t p = 0 ; p < N ; p++){
                acc[i*N + p] += phase[i*N + p];
        }
    }
    
    // The pointer to limb_0(b)
    PolyBiv* b_0 = res_ct + k*N;

    // For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[s_j * limb_i(a_j)]
    vec_znx_normalize_base2k_p(module, kappa, b_0, l, N*(k+1), acc, l, N);
    
    free(acc);

    return 0;
}

int ggsw_secret_encrypt(GGSWCtParams* enc_params,
                        GGSWCiphertext* res,           
                        GGSWSecretKeyDFT* sk_dft,             
                        PolyUniv* msg_univ)
{
    // GGSW & GLWE parameters
    GGSWCtParams* params_ggsw = enc_params;
    GLWECtParams* params_glwe = enc_params->params;
    uint64_t N = params_glwe->N;
    
    MODULE* module = new_module_info(N,FFT64);
    
    // Computes DFT(msg)
    PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe));
    if(msg_univ_dft == NULL){
        delete_module_info(module);
        perror("Malloc failed.");
        return -1;
    }
    vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, N);
    
    // We'll store DFT(msg * sk_j)
    PolyUnivDFT* phase_univ_dft = malloc(poly_univ_bytes(params_glwe));
    if (phase_univ_dft == NULL){
        delete_module_info(module);
        free(msg_univ_dft);  
        return -1;
    }

    // We'll store -msg * sk_j
    PolyUniv* phase_univ = malloc(poly_univ_bytes(params_glwe));
    if(phase_univ ==  NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft);
        perror("Malloc failed.");
        return -1;
    }

    // We'll store -msg * sk_j / (2^kappa_tilde)^i
    double* phase_univ_Rnx = malloc(poly_univ_bytes(params_glwe));
    if(phase_univ_Rnx ==  NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft); free(phase_univ);
        perror("Malloc failed.");
        return -1;
    }

    // We'll store the base-2^kappa decomposition of the phase = -m * sk_j / (2^kappa_tilde)^i 
    // and return in DFT space
    PolyBiv* phase_biv = malloc(poly_biv_bytes(params_glwe));
    if(phase_biv == NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft); free(phase_univ); free(phase_univ_Rnx);
        perror("Malloc failed.");
        return -1;
    }

    for (int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++){
        for (int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++){

            // The pointer to bivGLWE(-m * sk_j / (2^kappa_tilde)^i)
            VecBiv* ct_biv = ggsw_Sj_Yti(res, i, j);
            
            // Computes DFT(msg * sk_j)
            mult_vec_znx_dft(module, phase_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);
            
            // Computes -DFT(msg * sk_j)
            for(int64_t p = 0 ; p < N ; p++){
                phase_univ_dft[p] = -1 * phase_univ_dft[p];  
            }

            // Computes -msg * sk_j
            vec_znx_idft_p(module, phase_univ, 1, phase_univ_dft, 1);

            // Computes -msg * sk_j / (2^kappa_tilde)^i 
            for(int64_t p = 0 ; p < N ; p++){
                phase_univ_Rnx[p] = ldexp((double)phase_univ[p], -params_ggsw->kappa_tilde * i); 
            }

            // Compute the base-2^kappa decomposition of -m * sk_j / (2^kappa_tilde)^i 
            // and return in DFT space
            univ_to_biv(params_glwe, phase_biv, phase_univ_Rnx);

            // Computes the phase -m * sk_j / (2^kappa_tilde)^i + err
            add_error(params_glwe, phase_biv, phase_biv);

            #ifdef WITH_Y0 
            if (glwe_secret_masking_ggsw_lib(module, res->params->params, 
                                 ct_biv, sk_dft, phase_biv) < 0){
                return -1;
            }
            
            #endif 
            #ifndef WITH_Y0
            if (glwe_secret_masking(res->params->params, ct_m_sj_Yi,
                                 module, sk_dft, sk->size, 
                                 m_sk_j_dft) < 0){
                return -1;
            }
            #endif 
        }
    }

    delete_module_info(module);
    free(msg_univ_dft); free(phase_univ_dft); free(phase_univ); free(phase_univ_Rnx);

    return 0;
}

void ggsw_external_product(GLWECiphertext* res,  // result
                           GLWECiphertext* ct_glwe,  // GLWE ciphertext
                           GGSWCiphertext* ct_ggsw   // GGSW ciphertext
){
    uint64_t N = res->params->N;
    MODULE* module = new_module_info_p(N);

    // The bivGGSW ciphertext ct_ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs 
    // The bivGLWE ciphertext ct_glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
    // As the result of the vector-matrix product ct_glwe * ct_ggsw, 
    // the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
    uint64_t nrows = ct_ggsw->params->n_limbs_tilde;
    uint64_t ncols = ct_ggsw->params->params->n_limbs;

    // Computes the GGSW ciphertext in DFT space
    MatBivDFT* pmat = malloc(ggsw_coef_number(ct_ggsw->params) * sizeof(double));
    vec_znx_dft_p(module, pmat, nrows*ncols, ct_ggsw->mat, nrows*ncols, N);

    MatBivDFT* pvec = malloc(glwe_coef_number(res->params)*sizeof(double));

    // Computes ExternalProduct(ct_glwe, ct_ggsw)
    vmp_apply_dft_p(module, pvec, ncols, 
                            ct_glwe->vec, nrows, N,
                            pmat, nrows, ncols);

    vec_znx_idft_p(module, res->vec, ncols, pvec, ncols);

    free(pmat); free(pvec);                        
    delete_module_info_p(module);
}

//! GGSW IN DFT PART (begin)

void add_error_dft(GLWECtParams* enc_params,
                  PolyBivDFT* res_dft,
                  PolyBivDFT* phase_dft)
{
    MODULE* module = new_module_info(enc_params->N, FFT64);
    
    // Compute a random error in DFT space
    PolyBivDFT* err_dft = new_normal_random_biv_poly_dft(module, enc_params);
    
    // Add the error in DFT space
    add_biv_poly_dft(enc_params, phase_dft, enc_params->N, phase_dft, enc_params->N, err_dft, enc_params->N);
}

/**
 * @brief Computes Sum_j{0,k-1}[s_j * a_j].
 * 
 * @param module 
 * @param params 
 * @param res 
 * @param as 
 * @param sk_dft 
 * @return int 
 */
// TODO
int add_mult_dft_ggsw(const MODULE* module, GLWECtParams* params,
                 PolyBiv* res, VecBiv* as, GGSWSecretKeyDFT* sk_dft
){
    uint64_t N = params->N;
    uint64_t k = params->k;
    uint64_t l = poly_biv_size(params);

    // Computes Sum_j{0,k-1}[s_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes DFT(s_j) * DFT(a_j)
        PolyBivDFT* as_j_dft = calloc(poly_biv_coef_number_dft(params), 2*sizeof(double)); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, as + j*N, l, (k+1)*N); 
        
        // Computes s_j * a_j
        PolyBiv* as_j = calloc(poly_biv_coef_number(params), sizeof(int64_t)); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // And adds it to acc_j
        for(int64_t p = 0 ; p < N*l ; p++){
            res[p] += as_j[p];
        }
        free(as_j_dft); free(as_j);
    }
}

int glwe_secret_masking_dft(GLWECtParams* enc_params, 
                            const MODULE* module, 
                            VecBivDFT* res_dft,
                            GGSWSecretKeyDFT* sk_dft, 
                            PolyBivDFT* phase_dft)
{
    uint64_t N = enc_params->N;
    uint64_t k = enc_params->k;
    uint64_t kappa = enc_params->kappa;
    uint64_t n_limbs = enc_params->n_limbs;
    uint64_t l = n_limbs / (k+1);

    // Temporary bivGLWE ciphertext 
    VecBiv* tmp_ct = malloc(N*l*(k+1)*sizeof(int64_t));
    if(tmp_ct == NULL){
        perror("Malloc failed.");
        return -1;
    }
    
    // TODO coeff between 0 and (2^kappa)
    if (inplace_uniform_random_vec(k * N, tmp_ct, l, (k + 1) * N, kappa) < 0 )
    {
        free(tmp_ct);
        return -1;
    }
    
    // acc_(j+1) = acc_j + (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    
    // Computes Sum_j{0,k-1}[resVec_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of the secret key sk_dft
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        
        // Computes resVec_j_dft = (DFT(s_j) * limb_1(a_j) , ... , DFT(s_j) * limb_l(a_j))
        // TODO : can I only use one resVec_j, defined before the loop?
        PolyBivDFT* resVec_j_dft = new_vec_znx_dft_p(module, l); 
        svp_apply_dft_p(module, resVec_j_dft, l, sk_j_univ_dft, tmp_ct + j*N, l, (k+1)*N); 
        
        // Computes resVec_j in Zn[X,Y] space
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
    
    // The pointer to limb_0(b) in Zn[X,Y]
    PolyUniv* b_0_univ = tmp_ct + k*N;

    // For each i in {0,l} limb_i(b) = acc_i = Sum_j{0,k-1}[s_j * limb_i(a_j)]
    // Then b is normalized
    vec_znx_normalize_base2k_p(module, kappa, b_0_univ, l, N*(k+1), acc, l, N);
    
    // Computes tmp_ct in DFT space
    vec_znx_dft_p(module, res_dft, l*(k+1), tmp_ct, l*(k+1), N);
    
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

int ggsw_secret_encrypt_dft(GGSWCtParams* enc_params,
                            GGSWCiphertextDFT* res_dft,
                            GGSWSecretKeyDFT* sk_dft,
                            PolyUniv* msg_univ)
{
    // GGSW & GLWE parameters
    GGSWCtParams* params_ggsw = enc_params;
    GLWECtParams* params_glwe = enc_params->params;
    uint64_t N = params_glwe->N;

    // Prepare sk and m
    MODULE* module = new_module_info(N,FFT64);
    
    // Computes DFT(msg)
    PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe));
    if(msg_univ_dft == NULL){
        perror("Malloc failed.");
        delete_module_info(module);
        return -1;
    }
    vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, N);

    // Will store DFT(msg * sk_j)
    PolyUnivDFT* phase_univ_dft = malloc(poly_univ_bytes(params_glwe));
    if (phase_univ_dft == NULL){
        delete_module_info(module);
        free(msg_univ_dft);  
        return -1;
    }

    // Will store -msg * sk_j
    PolyUniv* phase_univ_inZ = malloc(poly_univ_bytes(params_glwe));
    if(phase_univ_inZ ==  NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft);
        perror("Malloc failed.");
        return -1;
    }

    // Will store -msg * sk_j / (2^kappa_tilde)^i
    double* phase_univ_inR = malloc(poly_univ_bytes(params_glwe));
    if(phase_univ_inR ==  NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft); free(phase_univ_inZ);
        perror("Malloc failed.");
        return -1;
    }

    // We'll store the base-2^kappa decomposition of the phase = -m * sk_j / (2^kappa_tilde)^i 
    PolyBiv* phase = malloc(poly_biv_bytes(params_glwe));
    if(phase == NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft); free(phase_univ_inZ); free(phase_univ_inR);
        perror("Malloc failed.");
        return -1;
    }

    // We'll store DFT(-m * sk_j / (2^kappa_tilde)^i)
    PolyBivDFT* phase_dft = malloc(poly_biv_bytes(params_glwe));
    if(phase_dft == NULL){
        delete_module_info(module);
        free(msg_univ_dft); free(phase_univ_dft); free(phase_univ_inZ); free(phase_univ_inR); free(phase);
        perror("Malloc failed.");
        return -1;
    }

    for (int64_t i = 0 ; i < nb_partials(params_ggsw) ; i++)
    {
        for (int64_t j = 0 ; j < nb_rows_per_partial(params_ggsw) ; j++)
        {
            // The pointer to bivGLWE(-m * s_j * Y^i) in DFT space
            VecBivDFT* ct_biv_dft = ggsw_Sj_Yti_dft(res_dft, i, j);
            
            // Computes DFT(-m * sk_j)
            mult_vec_znx_dft(module, phase_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);

            // Computes -msg * sk_j
            vec_znx_idft_p(module, phase_univ_inZ, 1, phase_univ_dft, 1);

            // Computes -msg * sk_j / (2^kappa_tilde)^i 
            for(int64_t p = 0 ; p < N ; p++){
                phase_univ_inR[p] = ldexp((double)phase_univ_inZ[p], -params_ggsw->kappa_tilde * i); 
            }

            // Compute the base-2^kappa decomposition of -m * sk_j / (2^kappa_tilde)^i 
            univ_to_biv(params_glwe, phase, phase_univ_inR);

            // Computes the phase -m * sk_j / (2^kappa_tilde)^i + err
            add_error(params_glwe, phase, phase);

            // Computes DFT(-m * sk_j / (2^kappa_tilde)^i + err)
            vec_znx_dft_p(module, phase_dft, poly_biv_size(params_glwe), phase, poly_biv_size(params_glwe), N);

            #ifdef WITH_Y0 
            // Computes bivGLWE(-m * s_j / (2^kappa_tilde)^i)
            if (glwe_secret_masking_dft(res_dft->params->params, module, ct_biv_dft,
                                  sk_dft, phase_dft) < 0)
            {
                delete_module_info(module);
                free(msg_univ_dft); free(phase_univ_dft); free(phase_univ_inZ); free(phase_univ_inR); free(phase); free(phase_dft);
                return -1;
            }

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
    delete_module_info(module);
    free(msg_univ_dft); free(phase_univ_dft); free(phase_univ_inZ); free(phase_univ_inR); free(phase); free(phase_dft);

    return 0;
}

void ggsw_external_product_dft(GLWECiphertextDFT* res_dft,  // result
                               GLWECiphertextDFT* ct_glwe_dft,  // GLWE ciphertext
                               GGSWCiphertextDFT* ct_ggsw_dft   // GGSW ciphertext
){
    uint64_t N = res_dft->params->N;
    MODULE* module = new_module_info_p(N);

    // The bivGGSW ciphertext ct_ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs 
    // The bivGLWE ciphertext ct_glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
    // As the result of the vector-matrix product ct_glwe * ct_ggsw, 
    // the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
    uint64_t nrows = ct_ggsw_dft->params->n_limbs_tilde;
    uint64_t ncols = ct_ggsw_dft->params->params->n_limbs;

    // Computes ExternalProduct(ct_glwe, ct_ggsw)
    vmp_apply_dft_to_dft_p(module, res_dft->pvec, ncols, 
                                   ct_glwe_dft->pvec, nrows, 
                                   ct_ggsw_dft->pmat, nrows, ncols);
}