#include "core/ggsw/ggsw.h"
#include "rng.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

#define NBASE 32
#define KBASE 1
#define KAPPABASE 4
#define NLIMBSBASE (KBASE + 1)*1
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -(LBASE/2 + 1)*KAPPABASE 

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE (K_TILDEBASE + 1)*1
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)
#define SIGMA_TILDEBASE -3


void printf_glwe(VecBiv* ct_glwe, GLWECtParams* params_glwe){
    printf("\n");
    // GLWE parameters
    uint64_t N = params_glwe->N;
    uint64_t k = params_glwe->k;
    uint64_t l = params_glwe->n_limbs / (params_glwe->k + 1);
    for(int64_t j = 0 ; j < params_glwe->k ; j++){
        printf("\na_%ld", j);
        printf_poly_biv(ct_glwe + j*N, (k+1)*N, N, l);
        printf("\n");
    }
    printf("\nb");
    printf_poly_biv(ct_glwe + k*N, (k+1)*N, N, l);
    printf("\n");
}

void printf_ggsw(MatBiv* ct_ggsw, GGSWCtParams* params_ggsw){
    printf("\n");
    // GGSW parameters
    uint64_t k_tilde = params_ggsw->k_tilde;
    uint64_t l_tilde = params_ggsw->n_limbs_tilde / (params_ggsw->k_tilde + 1);

    // GLWE parameters
    uint64_t N = params_ggsw->params_glwe->N;
    uint64_t k = params_ggsw->params_glwe->k;
    uint64_t l = params_ggsw->params_glwe->n_limbs / (params_ggsw->params_glwe->k + 1);

    for(int64_t i = 1 ; i <= l_tilde ; i++){
        for(int64_t j = 0 ; j < k_tilde; j++){
            printf("\nBivGLWE(-m * sk_%ld / 2^(kappa_tilde*%ld))", j, i+1);
            printf_glwe(ct_ggsw + (i-1)*(k_tilde + 1)*(k+1)*N*l + j*(k+1)*N*l, params_ggsw->params_glwe);
        }
        printf("\nBivGLWE(m / 2^(kappa_tilde*%ld))", i);
        printf_glwe(ct_ggsw + (i-1)*(k_tilde + 1)*(k+1)*N*l + k_tilde*(k+1)*N*l, params_ggsw->params_glwe);
    }
}

Test(smp_prepare, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    PolyBiv* pol = new_uniform_random_vec(NBASE, 3);

    SVP_PPOL* ppol = malloc(NBASE*sizeof(double));
    svp_prepare_p(module, (double*)ppol, pol);

    PolyBivDFT* pol_dft = malloc(NBASE*sizeof(double));
    vec_znx_dft_p(module, pol_dft, 1, pol, 1, NBASE);

    printf("ppol\n");
    printf_poly_univ_RnX((double*)ppol, NBASE);
    printf("\n");

    printf("pol_dft\n");
    printf_poly_univ_RnX(pol_dft, NBASE);
    printf("\n");

    PolyBiv* ppol_out = malloc(NBASE * sizeof(double));
    vec_znx_idft_p(module, ppol_out, 1, (double*)ppol, 1);

    PolyBiv* pol_dft_out = malloc(NBASE * sizeof(double));
    vec_znx_idft_p(module, pol_dft_out, 1, pol_dft, 1);

    printf("ppol_out\n");
    printf_poly_univ_ZnX(ppol_out, NBASE);
    printf("\n");

    printf("pol_dft_out\n");
    printf_poly_univ_ZnX(pol_dft_out, NBASE);
    printf("\n");

    delete_module_info(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

Test(vmp_prepare, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    uint64_t nrows = params->n_limbs_tilde;
    uint64_t ncols = params->params_glwe->n_limbs;

    MatBiv* mat = new_uniform_random_vec(ggsw_coef_number(params), 3);

    VMP_PMAT* pmat = malloc(ggsw_bytes(params));
    vmp_prepare_contiguous_p(module, (double*)pmat, mat, nrows, ncols);

    MatBivDFT* mat_dft = malloc(ggsw_bytes(params));
    vec_znx_dft_p(module, mat_dft, ggsw_size(params), mat, ggsw_size(params), NBASE);

    // printf("pmat\n");
    // for(int64_t i = 0 ; i < ncols ; i++){
    //     for(int64_t j = 0 ; j < nrows ; j++){
    //         printf("\n i = %ld j = %ld \n", i, j);
    //         printf_poly_univ_RnX((double*)pmat + i*nrows*NBASE + j*NBASE, NBASE);
    //         printf("\n");
    //     }
    // }

    // printf("mat_dft\n");
    // for(int64_t i = 0 ; i < nrows ; i++){
    //     for(int64_t j = 0 ; j < ncols ; j++){
    //         printf("\n i = %ld j = %ld \n", i, j);
    //         printf_poly_univ_RnX(mat_dft + i*ncols*NBASE + j*NBASE, NBASE);
    //         printf("\n");
    //     }
    // }

    MatBiv* pmat_out = malloc(ggsw_bytes(params));
    vec_znx_idft_p(module, pmat_out, ggsw_size(params), (double*)pmat, ggsw_size(params));

    MatBiv* mat_dft_out = malloc(ggsw_bytes(params));
    vec_znx_idft_p(module, mat_dft_out, ggsw_size(params), mat_dft, ggsw_size(params));

    // printf("pmat_out\n");
    // for(int64_t i = 0 ; i < nrows ; i++){
    //     for(int64_t j = 0 ; j < ncols ; j++){
    //         printf("\n i = %ld j = %ld \n", i, j);
    //         printf_poly_univ_ZnX(pmat_out + i*ncols*NBASE + j*NBASE, NBASE);
    //         printf("\n");
    //     }
    // }

    // printf("mat_dft_out\n");
    // for(int64_t i = 0 ; i < ncols ; i++){
    //     for(int64_t j = 0 ; j < nrows ; j++){
    //         printf("\n i = %ld j = %ld \n", i, j);
    //         printf_poly_univ_ZnX(mat_dft_out + i*nrows*NBASE + j*NBASE, NBASE);
    //         printf("\n");
    //     }
    // }

    delete_module_info(module);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

/**
 * @brief Tests ggsw_secret_encrpyt
 * 
 */
Test(ggsw_secret_encrypt, works)
{
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = ldexp(1.0, -(LBASE/2 + 1)*KAPPABASE);

    // The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99% chance
    double err_length = ldexp(1.0,-(LBASE/2)*KAPPABASE) + 3*sigma;
    cr_log_info("error length = %e", err_length);

    GLWECtParams* params_glwe_for_ggsw = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe_for_ggsw, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info_p(NBASE);
    
    GGSWCiphertext* ct_ggsw = new_ggsw(params_ggsw, NULL);
    GGSWSecretKeyDFT* sk_dft = new_uniform_ggsw_secret_key_dft(module, KBASE, 3);

    // Message uniformly drawn
    PolyUniv* msg_univ = new_uniform_random_vec(NBASE, KAPPABASE*LBASE);
    PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
    vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, NBASE);

    // Computes a bivGGSW(msg)
    ggsw_secret_encrypt(module, params_ggsw, ct_ggsw, sk_dft, msg_univ);

    for(int64_t i = 1 ; i <= L_TILDEBASE ; i++){
        for(int64_t j = 0 ; j < K_TILDEBASE ; j++){
            // The pointer to bivGLWE(-m * sk_j / (2^kappa_tilde)^(i+1))  
            VecBiv* ct_glwe = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            memcpy(ct_glwe, ct_ggsw->mat + ((i-1)*(K_TILDEBASE+1) + j)*glwe_coef_number(params_glwe_for_ggsw), glwe_bytes(params_glwe_for_ggsw));

            // Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1)) + err
            PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            glwe_secret_demasking_ggsw_lib(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe);

            // Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1) + err in RnX
            PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);

            // Computes by hand the phase = -m * sk_j / 2^kappa_tilde*(i+1) 
            PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            PolyUnivDFT* m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
            PolyUniv* m_skj_univ = malloc(poly_univ_bytes(params_glwe_for_ggsw));
            
            // Computes DFT(msg * sk_j)
            mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);
            
            // Computes -msg * sk_j
            for(int64_t p = 0 ; p < NBASE ; p++){
                m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];  
            } 
            vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1);
           
            // Computes -msg * sk_j / (2^kappa_tilde)^(i+1) 
            for(int64_t p = 0 ; p < NBASE ; p++){
                phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -(params_ggsw->kappa_tilde*(i+1))); 
            }

            // Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase, 
            // are only different by an error of approximation and a gaussian error
            for(int64_t p = 0 ; p < NBASE ; p++){
                double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
                double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] + floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);
                double err_length = ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE);

                int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

                cr_assert(cond,
                    "Equality failed at p = %ld with : \n-(msg * sk_j)[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \n-(msg * sk_j)[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf", 
                p, p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p], ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE), 
                   p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]), ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE));
            }

            free(phase_univ_RnX_computed); free(m_skj_univ_dft); free(m_skj_univ);
            free(ct_glwe); free(phase); free(phase_univ_RnX);
        }
        // The pointer to bivGLWE(msg / (2^kappa_tilde)^(i+1))  
        VecBiv* ct_glwe = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
        memcpy(ct_glwe, ct_ggsw->mat + ((i-1)*(K_TILDEBASE+1) + K_TILDEBASE)*glwe_coef_number(params_glwe_for_ggsw), glwe_bytes(params_glwe_for_ggsw));
        
        // Computes the phase = m/2^kappa_tilde + err
        PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
        glwe_secret_demasking_ggsw_lib(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe);
    
        // Computes the phase in Rn[X]
        PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
        biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);
    
        // Computes by hand the phase = msg / (2^kappa_tilde)^(i+1) 
        PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
    
        // Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase, 
        // are only different by an error of approximation and a gaussian error
        for(int64_t p = 0 ; p < NBASE ; p++){
            double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
            double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] + floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);
            
            int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

            cr_assert(cond, "Equality failed at p = %ld with : \nmsg[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \nmsg[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf", 
            p, p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p], ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE), 
               p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]), ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE));
        
        }
        free(ct_glwe); free(phase); free(phase_univ_RnX); free(phase_univ_RnX_computed);
    }

    free(msg_univ); free(msg_univ_dft);
    delete_module_info_p(module);
    delete_ggsw(ct_ggsw);
    delete_ggsw_secret_key_dft(sk_dft);
    delete_glwe_ct_params(params_glwe_for_ggsw);
    delete_ggsw_ct_params(params_ggsw);
}

/**
 * @brief Construct a new Test object
 * 
 */
Test(ggsw_secret_encrypt_dft, works)
{
    // GLWE and GGSW parameters. This set of GLWE parameters is for GGSW ciphertext
    double sigma = ldexp(1.0, -(LBASE/2 + 1)*KAPPABASE);

    // The decryption of a bivGLWE(m) should give m_dec = m + err, and |m_dec - m| <= 2^(-kappa*l) + 3*sigma with a 99% chance
    double err_length = ldexp(1.0,-(LBASE/2)*KAPPABASE) + 3*sigma;
    cr_log_info("error length = %e", err_length);

    GLWECtParams* params_glwe_for_ggsw = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, sigma);
    GGSWCtParams* params_ggsw = new_ggsw_ct_params(params_glwe_for_ggsw, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info_p(NBASE);
    
    GGSWCiphertextDFT* ct_ggsw_dft = new_ggsw_dft(params_ggsw, NULL);
    GGSWSecretKeyDFT* sk_dft = new_uniform_ggsw_secret_key_dft(module, KBASE, 3);

    // Message uniformly drawn
    PolyUniv* msg_univ = new_uniform_random_vec(NBASE, KAPPABASE*LBASE);
    PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
    vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, NBASE);

    // Computes a bivGGSW(msg)
    ggsw_secret_encrypt_dft(module, params_ggsw, ct_ggsw_dft, sk_dft, msg_univ);

    for(int64_t i = 1 ; i <= L_TILDEBASE ; i++){
        for(int64_t j = 0 ; j < K_TILDEBASE ; j++){
            // The pointer to DFT(bivGLWE(-m * sk_j / (2^kappa_tilde)^(i+1)))
            VecBivDFT* ct_glwe_dft = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            memcpy(ct_glwe_dft, ct_ggsw_dft->mat + ((i-1)*(K_TILDEBASE+1) + j)*glwe_coef_number(params_glwe_for_ggsw), glwe_bytes(params_glwe_for_ggsw));

            // Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1)) + err
            PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            glwe_secret_demasking_ggsw_lib_dft(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe_dft);

            // Computes the phase = -m * sk_j / (2^kappa_tilde)^(i+1) + err in RnX
            PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);

            // Computes by hand the phase = -m * sk_j / 2^kappa_tilde*(i+1) 
            PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
            PolyUnivDFT* m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe_for_ggsw));
            PolyUniv* m_skj_univ = malloc(poly_univ_bytes(params_glwe_for_ggsw));
            
            // Computes DFT(msg * sk_j)
            mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);
            
            // Computes -msg * sk_j
            for(int64_t p = 0 ; p < NBASE ; p++){
                m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];  
            } 
            vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1);
           
            // Computes -msg * sk_j / (2^kappa_tilde)^(i+1) 
            for(int64_t p = 0 ; p < NBASE ; p++){
                phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -(params_ggsw->kappa_tilde*(i+1))); 
            }

            // Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase, 
            // are only different by an error of approximation and a gaussian error
            for(int64_t p = 0 ; p < NBASE ; p++){
                double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
                double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] + floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);
                double err_length = ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE);

                int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

                cr_assert(cond,
                    "Equality failed at p = %ld with : \n-(msg * sk_j)[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \n-(msg * sk_j)[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf", 
                p, p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p], ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE), 
                   p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]), ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE));
            }

            free(phase_univ_RnX_computed); free(m_skj_univ_dft); free(m_skj_univ);
            free(ct_glwe_dft); free(phase); free(phase_univ_RnX);
        }
        // The pointer to bivGLWE(msg / (2^kappa_tilde)^(i+1))  
        VecBivDFT* ct_glwe_dft = calloc(glwe_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
        memcpy(ct_glwe_dft, ct_ggsw_dft->mat + ((i-1)*(K_TILDEBASE+1) + K_TILDEBASE)*glwe_coef_number(params_glwe_for_ggsw), glwe_bytes(params_glwe_for_ggsw));
        
        // Computes the phase = m/2^kappa_tilde + err
        PolyBiv* phase = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
        glwe_secret_demasking_ggsw_lib_dft(module, params_glwe_for_ggsw, phase, sk_dft, ct_glwe_dft);
    
        // Computes the phase in Rn[X]
        PolyUnivRnX* phase_univ_RnX_computed = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
        biv_to_univ(params_glwe_for_ggsw, phase_univ_RnX_computed, phase);
    
        // Computes by hand the phase = msg / (2^kappa_tilde)^(i+1) 
        PolyUnivRnX* phase_univ_RnX = calloc(poly_biv_coef_number(params_glwe_for_ggsw), sizeof(int64_t));
    
        // Assures that the difference between the phase = msg / (2^kappa_tilde)^(i+1) and the computed phase, 
        // are only different by an error of approximation and a gaussian error
        for(int64_t p = 0 ; p < NBASE ; p++){
            double diff_1 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p];
            double diff_2 = phase_univ_RnX[p] - floor(phase_univ_RnX[p]) - phase_univ_RnX_computed[p] + floor(phase_univ_RnX_computed[p]) + ceil(phase_univ_RnX_computed[p]);
            
            int cond = (diff_1 <= err_length || diff_1 >= -err_length) || (diff_2 <= err_length || diff_2 >= -err_length);

            cr_assert(cond, "Equality failed at p = %ld with : \nmsg[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf, \nmsg[%ld] / (2^kappa_tilde)^%ld = %lf and phase_univ_RnX_computed[%ld] = %lf and error_length = %lf", 
            p, p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p], ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE), 
               p, K_TILDEBASE, i+1, phase_univ_RnX[p] - floor(phase_univ_RnX[p]), p, phase_univ_RnX_computed[p] - floor(phase_univ_RnX_computed[p]) - ceil(phase_univ_RnX_computed[p]), ldexp(1.0,-(LBASE/2)*KAPPABASE) + ldexp(1.0, -LBASE*KAPPABASE));
        
        }
        free(ct_glwe_dft); free(phase); free(phase_univ_RnX); free(phase_univ_RnX_computed);
    }

    free(msg_univ); free(msg_univ_dft);
    delete_module_info_p(module);
    delete_ggsw_dft(ct_ggsw_dft);
    delete_ggsw_secret_key_dft(sk_dft);
    delete_glwe_ct_params(params_glwe_for_ggsw);
    delete_ggsw_ct_params(params_ggsw);
}