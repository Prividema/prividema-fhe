#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ct_params.h"
#include "spqlios_alias.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>

#define NBASE 4
#define KBASE 8
#define KAPPABASE 32
#define NLIMBSBASE 180
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

Test(vec_znx_normalize_base2k_p, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    MODULE* module = new_module_info_p(params->N);
    
    PolyBiv* a = malloc(LBASE * NBASE * sizeof(int64_t));
    PolyBiv* res = malloc(LBASE * NBASE * sizeof(int64_t));
    vec_znx_normalize_base2k_p(module, KAPPABASE, res, LBASE, NBASE, a, LBASE, NBASE);

    free(a); free(res);
    delete_glwe_ct_params(params);
    delete_module_info_p(module);
}

Test(vmp_product, basic){
    MODULE* module = new_module_info_p(NBASE);
    uint64_t nrows = 2 ;
    uint64_t ncols = 1 ;
    double* pmat = new_vmp_pmat_p(module, nrows, ncols);

    // mat = (1 + X)
    //       (1 + 2X)
    MatUniv* mat = calloc(nrows * ncols * NBASE , sizeof(int64_t));

    // P = 1 + X
    mat[0] = 1;
    mat[1] = 1;

    // Q = 1 + 2X
    mat[NBASE + 0] = 1;
    mat[NBASE + 1] = 2;

    vmp_prepare_contiguous_p(module, pmat, mat, nrows, ncols);

    // Vec = (4 , 1 + 8X)
    VecUniv* a = calloc(nrows * NBASE, sizeof(int64_t));

    a[0] = 4;
    
    a[NBASE + 0] = 1;
    a[NBASE + 1] = 8;

    // Result 
    VecUnivDFT* res_dft = calloc(ncols * NBASE, sizeof(double));
    vmp_apply_dft_p(module, res_dft, 1, a, 2, NBASE, pmat, nrows, ncols);

    VecUniv* res = calloc(ncols * NBASE, sizeof(int64_t));
    vec_znx_idft_p(module, res, ncols, res_dft, ncols);

    for(int64_t i = 0 ; i < ncols ; i++){
        for(int64_t p = 0 ; p < NBASE ; p++)
            cr_log_info("%ld X^%ld", res[i*NBASE + p], p);
    }

    free(a); free(mat); free(res_dft); free(res); free(pmat);
    delete_module_info_p(module);
}
