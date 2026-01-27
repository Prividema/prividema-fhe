#include "bivariate_polynomial.h"
#include "core/glwe/glwe.h"
#include "common/distributions.h"
#include "common/spqlios_alias.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>
#include <stdio.h>

#define NBASE 8
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -7

/**
 * @brief Test univ_to_biv
 */
Test(univ_to_biv, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, ldexp(1.0, SIGMABASE));
    MODULE* module = new_module_info_p(NBASE);

    double* pol_univ = malloc(poly_univ_bytes(params));
    new_normal_random_vec(NBASE, pol_univ, 1, NBASE, 0.0, params->sigma);

    PolyBiv* pol_biv = malloc(poly_biv_bytes(params));
    univ_to_biv(params, pol_biv, pol_univ);
    
    for(int64_t p = 0 ; p < NBASE ; p++){
        printf("\n\nA_p : %e X^%ld \n", pol_univ[p], p);
        for(int64_t i = 1 ; i < LBASE ; i++)
            printf(" %ld Y^%ld ", pol_biv[i*NBASE + p], i);
    }

    for(int64_t p = 0 ; p < NBASE ; p++){
        double acc = 0; 
        for(int64_t i = 1 ; i < LBASE ; i++){
            acc += ldexp((double)pol_biv[i * NBASE + p], -i * KAPPABASE);
        }
        cr_log_info("acc %lf pol %lf p %ld", acc, pol_univ[p], p);
        cr_assert(epsilon_eq(dbl, acc - floor(acc) - pol_univ[p] + floor(pol_univ[p]), 0, ldexp(1.0,-(LBASE-1)*KAPPABASE))
                  , "acc %lf pol %lf p %ld", acc - floor(acc), - pol_univ[p] + floor(pol_univ[p]), p);
    }
    free(pol_univ); free(pol_biv);
    delete_glwe_ct_params(params);
    delete_module_info_p(module);
}
