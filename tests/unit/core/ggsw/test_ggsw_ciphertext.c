#include "ggsw_ciphertext.h"
#include "rng.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 4
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE (KBASE + 1)*2
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

#define K_TILDEBASE 1
#define KAPPA_TILDEBASE 4
#define NLIMBS_TILDEBASE (K_TILDEBASE+1)*1
#define L_TILDEBASE NLIMBS_TILDEBASE/(K_TILDEBASE+1)

//! COMMON PART (begin)

// Test ggsw_size
Test(ggsw_size, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, ggsw_size(params), NLIMBS_TILDEBASE * NLIMBSBASE));

    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test ggsw_bytes
Test(ggsw_bytes, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, ggsw_bytes(params), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE * 8));

    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}


//! GGSW Part (begin)

// Test ggsw_coef_number
Test(ggsw_coef_number, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    cr_assert(eq(i64, ggsw_coef_number(params), NLIMBS_TILDEBASE * NLIMBSBASE * NBASE));
    
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test new_ggsw
Test(new_ggsw, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct = new_ggsw(params, NULL);

    cr_assert(eq(int, ct->mat != NULL, 1));
    cr_assert(eq(int, ct->params != NULL, 1));
    
    delete_ggsw(ct);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test ggsw_Sj_Yti
Test(ggsw_Sj_Yti, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct = new_ggsw(params, NULL);

    for(int64_t i = 1 ; i < nb_partials(params) ; i++)
        for(int64_t j = 0 ; j < K_TILDEBASE + 1 ; j++)
        {
            VecBiv* ct_mat_ij = ggsw_Sj_Yti(params, ct->mat, j, i);

            // Modify the two firsts coefficients of biGLWE(-m * sk_j / (2^kappa_tilde)^i).
            ct_mat_ij[0] = 1;
            ct_mat_ij[1] = 2;
            
            cr_assert(eq(i64, ct->mat[(i-1)*(K_TILDEBASE + 1)*NLIMBSBASE*NBASE + j*NLIMBSBASE*NBASE], 1));
            cr_assert(eq(i64, ct->mat[(i-1)*(K_TILDEBASE + 1)*NLIMBSBASE*NBASE + j*NLIMBSBASE*NBASE + 1], 2));
        }
    
    delete_ggsw(ct);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

// Test normalize_ggsw
Test(normalize_ggsw, basic){
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info(NBASE, FFT64);

    GGSWCiphertext* ct = new_ggsw(params, NULL);
    GGSWCiphertext* res = new_ggsw(params, NULL);

    normalize_ggsw(module, res, ct);

    delete_ggsw(ct);
    delete_ggsw(res);
    delete_glwe_ct_params(params_glwe);
    delete_ggsw_ct_params(params);
}

Test(add_ggsw, basic)
{
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);

    GGSWCiphertext* ct_l = new_ggsw(params, NULL);
    GGSWCiphertext* ct_r = new_ggsw(params, NULL);
    GGSWCiphertext* res = new_ggsw(params, NULL);

    // Draws uniformly the GGSW ciphertexts
    inplace_uniform_random_vec(NBASE, ct_l->mat, ggsw_size(params), NBASE, KAPPABASE-1);
    inplace_uniform_random_vec(NBASE, ct_r->mat, ggsw_size(params), NBASE, KAPPABASE-1);

    add_ggsw(res, ct_l, ct_r);
    
    uint64_t nb_rows = res->params->n_limbs_tilde;
    uint64_t nb_cols = res->params->params_glwe->n_limbs;
    uint64_t N = res->params->params_glwe->N;

    for (int64_t i = 0 ; i < nb_rows ; i++)
        for(int64_t j = 0 ; j < nb_cols ; j++)
            for (uint64_t p = 0 ; p < N ; p++)
                cr_assert(eq(res->mat[i*N*nb_cols + j*N + p], ct_l->mat[i*N*nb_cols + j*N + p] + ct_r->mat[i*N*nb_cols + j*N + p]));
    
    delete_ggsw(ct_l);
    delete_ggsw(ct_r);
    delete_ggsw(res);
    delete_ggsw_ct_params(params);
    delete_glwe_ct_params(params_glwe);
}

Test(const_mult_ggsw, without_normalization)
{
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info(NBASE, FFT64);
    
    PolyUniv* u = new_uniform_random_vec(NBASE, KAPPABASE-1);
    GGSWCiphertext* ct = new_ggsw(params, NULL);
    GGSWCiphertext* res = new_ggsw(params, NULL);
    
    // Draws uniformly the GGSW ciphertext 
    inplace_uniform_random_vec(NBASE, ct->mat, ggsw_size(params), NBASE, KAPPABASE-1);

    PolyUnivDFT* u_dft = malloc(NBASE * sizeof(double));
    vec_znx_dft_p(module, u_dft, 1, u, 1, NBASE);

    const_mult_ggsw(module, res, ct, u_dft, 0);

    for(int64_t ii = 1 ; ii <= nb_partials(params) ; ii++)
        for(int64_t jj = 0 ; jj < K_TILDEBASE + 1 ; jj++)
        {
            VecBiv* ct_mat_ii_jj = ggsw_Sj_Yti(params, ct->mat, jj, ii);
            VecBiv* res_mat_ii_jj = ggsw_Sj_Yti(params, res->mat, jj, ii);
            for(int64_t j = 0 ; j < KBASE + 1 ; j++)
                for(int64_t p = 0 ; p < NBASE ; p++)for(int64_t i = 1 ; i <= LBASE ; i++)
                {
                    for(int64_t i = LBASE ; i >= 1 ; i--)
                    {
                        PolyUniv* ct_ij = ct_mat_ii_jj + (i-1)*(KBASE+1)*NBASE + j*NBASE;
                        int64_t acc = 0;
                        for(uint64_t k = 0 ; k <= p; k++)
                        {
                            acc += u[k] * ct_ij[p-k];
                        } 
                        for(uint64_t k = p + 1; k < NBASE; k++)
                        {
                            acc += -u[k] * ct_ij[NBASE + p-k];
                        }  
                        cr_assert(eq(i64, res_mat_ii_jj[(i-1)*(KBASE+1)*NBASE + j*NBASE + p], acc));
                    }
                }
        }

    free(u); free(u_dft);
    delete_module_info(module);
    delete_ggsw(ct);
    delete_ggsw(res);
    delete_ggsw_ct_params(params);
    delete_glwe_ct_params(params_glwe);
}

Test(const_mult_ggsw, with_normalization)
{
    GLWECtParams* params_glwe = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GGSWCtParams* params = new_ggsw_ct_params(params_glwe, K_TILDEBASE, KAPPA_TILDEBASE, NLIMBS_TILDEBASE);
    MODULE* module = new_module_info(NBASE, FFT64);
    
    PolyUniv* u = new_uniform_random_vec(NBASE, KAPPABASE-1);
    GGSWCiphertext* ct = new_ggsw(params, NULL);
    GGSWCiphertext* res = new_ggsw(params, NULL);
    
    // Draws uniformly the GGSW ciphertext 
    inplace_uniform_random_vec(NBASE, ct->mat, ggsw_size(params), NBASE, KAPPABASE-1);

    PolyUnivDFT* u_dft = malloc(NBASE * sizeof(double));
    vec_znx_dft_p(module, u_dft, 1, u, 1, NBASE);

    const_mult_ggsw(module, res, ct, u_dft, 1);

    for(int64_t ii = 1 ; ii <= nb_partials(params) ; ii++)
        for(int64_t jj = 0 ; jj < K_TILDEBASE + 1 ; jj++)
        {
            VecBiv* ct_mat_ii_jj = ggsw_Sj_Yti(params, ct->mat, jj, ii);
            VecBiv* res_mat_ii_jj = ggsw_Sj_Yti(params, res->mat, jj, ii);
            
            for(int64_t j = 0 ; j < KBASE + 1 ; j++)
                for(int64_t p = 0 ; p < NBASE ; p++)for(int64_t i = 1 ; i <= LBASE ; i++)
                {
                    int64_t remainder = 0;
                    for(int64_t i = LBASE ; i >= 1 ; i--)
                    {
                        PolyUniv* ct_ij = ct_mat_ii_jj + (i-1)*(KBASE+1)*NBASE + j*NBASE;

                        int64_t acc = 0;
                        for(uint64_t k = 0 ; k <= p; k++)
                        {
                            acc += u[k] * ct_ij[p-k];
                        } 
                        for(uint64_t k = p + 1; k < NBASE; k++)
                        {
                            acc += -u[k] * ct_ij[NBASE + p-k];
                        }  

                        cr_assert(eq(i64, (res_mat_ii_jj[(i-1)*(KBASE+1)*NBASE + j*NBASE + p] - (acc + remainder)) % (1 << KAPPABASE), 0),
                                "Equality failed at j = %ld p = %ld i = %ld with acc = %ld reminder = %ld and res = %ld", 
                                j, p, i, acc, remainder, res_mat_ii_jj[(i-1)*(KBASE+1)*NBASE + j*NBASE + p]);
                        
                        remainder = acc >= 0 ? (acc + (1 << KAPPABASE-1)) / (1 << KAPPABASE) : (acc - (1 << KAPPABASE-1) + 1) / (1 << KAPPABASE);
                    }
                }
        }

    free(u); free(u_dft);
    delete_module_info(module);
    delete_ggsw(ct);
    delete_ggsw(res);
    delete_ggsw_ct_params(params);
    delete_glwe_ct_params(params_glwe);
}
// void const_mult_ggsw(GGSWCiphertext* res,  
//                      GGSWCiphertext* ct, 
//                      PolyUniv* u,
//                      int do_normalization)

// TODO