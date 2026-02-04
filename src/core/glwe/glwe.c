#include "glwe.h"
#include "rng.h"

//! GLWE PART (begin)
int add_mult(GLWECtParams* params, MODULE* module, 
             PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft
){
    // GLWE parameters
    uint64_t N = params->N;
    uint64_t k = params->k;
    uint64_t l = poly_biv_size(params);
    
    // Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of resp. the secret key and the bivGLWE ciphertext 
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        PolyBiv* a_j = ct + j*N;
        
        // Computes DFT(sk_j * a_j)
        PolyBivDFT* as_j_dft = malloc(poly_biv_bytes(params)); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k+1)*N); 
        
        // Computes sk_j * a_j
        PolyBiv* as_j = malloc(poly_biv_bytes(params)); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // Computes acc = acc - sk_j * a_j
        for(int64_t p = 0 ; p < N*l ; p++){
            res[p] += as_j[p];
        }
        free(as_j_dft);
        free(as_j);
    }
}

int glwe_secret_masking(GLWECiphertext* ct, 
                        GLWESecretKeyDFT* sk_dft,  
                        PolyBiv* phase  
){
    uint64_t N = ct->params->N;
    uint64_t k = ct->params->k;
    uint64_t kappa = ct->params->kappa;
    uint64_t l = poly_biv_size(ct->params);

    MODULE* module = new_module_info(N, FFT64);
    if (inplace_uniform_random_vec(k * N, ct->vec, l, (k + 1) * N, kappa) < 0) {
        delete_module_info_p(module);
        return -1;
    }
    
    // acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    if (!acc){
        perror("calloc failed");
        delete_module_info_p(module);
        return -1;
    }

    // Computes Sum_j{0,k-1}[s_j * a_j]
    if (add_mult(ct->params, module, acc, ct->vec, sk_dft) < 0){
        return -1;
    }

    // Add the phase to acc
    for(int64_t i = 0 ; i < l ; i++){    
        for(int64_t p = 0 ; p < N ; p++){
                acc[i*N + p] += phase[i*N + p];
        }
    }
    
    // The pointer to limb_0(b)
    PolyBiv* b_0 = ct->vec + k*N;

    // For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[s_j * limb_i(a_j)]
    vec_znx_normalize_base2k_p(module, kappa, b_0, l, N*(k+1), acc, l, N);
    
    free(acc);
    delete_module_info_p(module);

    return 0;
}

int sub_mult(GLWECtParams* params, MODULE* module, 
             PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft
){
    // GLWE parameters
    uint64_t N = params->N;
    uint64_t k = params->k;
    uint64_t l = poly_biv_size(params);
    
    // Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
    for(int64_t j = 0 ; j < k ; j++)
    {
        // The j-ème component of resp. the secret key and the bivGLWE ciphertext 
        PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j]; 
        PolyBiv* a_j = ct + j*N;
        
        // Computes DFT(sk_j * a_j)
        PolyBivDFT* as_j_dft = malloc(poly_biv_bytes(params)); 
        svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k+1)*N); 
        
        // Computes sk_j * a_j
        PolyBiv* as_j = malloc(poly_biv_bytes(params)); 
        vec_znx_idft_p(module, as_j, l, as_j_dft, l);

        // Computes acc = acc - sk_j * a_j
        for(int64_t p = 0 ; p < N*l ; p++){
            res[p] -= as_j[p];
        }
        free(as_j_dft);
        free(as_j);
    }
}

int glwe_secret_demasking(PolyBiv* res,  
                          GLWESecretKeyDFT* sk_dft, 
                          GLWECiphertext* ct 
){
    // GLWE parameters
    uint64_t N = ct->params->N;
    uint64_t k = ct->params->k;
    uint64_t l = poly_biv_size(ct->params);

    MODULE* module = new_module_info(N, FFT64);

    PolyBiv* acc = calloc(N*l, sizeof(int64_t)); 
    if (acc == NULL){
        perror("calloc failed");
        delete_module_info_p(module);
        return -1;
    }

    if (sub_mult(ct->params, module, acc, ct->vec, sk_dft) < 0){
        return -1;
    }

    // Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
    int64_t* b = ct->vec + N*k;
    add_biv_poly(ct->params, acc, N, b, N*(k+1), acc, N);
    
    // The phase in Zn[X,Y]
    vec_znx_normalize_base2k_p(module, ct->params->kappa, res, l, N, acc, l, N);

    free(acc);
    delete_module_info_p(module);
    
    return 0;
}

//! GLWE IN DFT SPACE PART (begin)

int glwe_secret_masking_dft(GLWECiphertextDFT* ct_dft,  
                            GLWESecretKeyDFT* sk_dft,   
                            PolyBivDFT* phase_dft   
){
    // GLWE parameters
    GLWECtParams* params = ct_dft->params;
    uint64_t k = params->k;
    uint64_t N = params->N;
    uint64_t kappa = params->kappa;
    uint64_t l = poly_biv_size(params);

    MODULE* module = new_module_info(N, FFT64);
    VecBiv* ct = calloc(glwe_coef_number(params), sizeof(int64_t));

    if (inplace_uniform_random_vec(k * N, ct, l, (k + 1) * N, kappa) < 0) {
        delete_module_info_p(module);
        return -1;
    }

    // acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
    PolyBiv* acc = calloc(N*l,sizeof(double)); 
    if (!acc){
        perror("calloc failed");
        delete_module_info_p(module);
        return -1;
    }

    // Computes Sum_j{0,k-1}[s_j * a_j]
    if (add_mult(params, module, acc, ct, sk_dft) < 0){
        return -1;
    }

    // The pointer to limb_0(b)
    PolyBiv* b_0 = ct + k*N;

    // For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[s_j * limb_i(a_j)]
    vec_znx_normalize_base2k_p(module, kappa, b_0, l, N*(k+1), acc, l, N);
    
    // Computes the GLWE ciphertext in DFT space
    vec_znx_dft_p(module, ct_dft->pvec, (k+1)*l, ct, (k+1)*l, N);

    // The phase 
    PolyBiv* phase = malloc(poly_biv_bytes(params));
    vec_znx_idft_p(module, phase, l, phase_dft, l);

    // Add the phase to the result ciphertext's b
    for(int64_t i = 0 ; i < l ; i++){    
        for(int64_t p = 0 ; p < N ; p++){
            ct_dft->pvec[i*N*(k+1) + N*k + p] += phase_dft[i*N + p];
        }
    }

    free(acc); free(ct); free(phase);
    delete_module_info_p(module);

    return 0;
}

int glwe_secret_demasking_dft(PolyBiv* res,  
                              GLWESecretKeyDFT* sk_dft, 
                              GLWECiphertextDFT* ct_dft 
){
    // GLWE parameters
    GLWECtParams* params = ct_dft->params;
    uint64_t N = params->N;
    uint64_t k = params->k;
    uint64_t l = poly_biv_size(params);

    MODULE* module = new_module_info(N, FFT64);
    
    // Computes the input ciphertext out of DFT space
    VecBiv* ct = calloc(glwe_coef_number(params), sizeof(int64_t));
    vec_znx_idft_p(module, ct, glwe_size(params), ct_dft->pvec, glwe_size(params));

    PolyBiv* acc = calloc(N*l, sizeof(int64_t)); 
    if (acc == NULL){
        perror("calloc failed");
        delete_module_info_p(module);
        return -1;
    }

    if (sub_mult(params, module, acc, ct, sk_dft) < 0){
        return -1;
    }

    // Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
    PolyBiv* b = ct + N*k;
    add_biv_poly(params, acc, N, b, N*(k+1), acc, N);
    
    // The phase in Zn[X,Y]
    vec_znx_normalize_base2k_p(module, params->kappa, res, l, N, acc, l, N);

    free(acc); free(ct);
    delete_module_info_p(module);
    
    return 0;
}
