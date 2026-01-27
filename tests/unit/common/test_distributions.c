#include "common/distributions.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include <stdio.h> 

#define NBASE 4
#define KBASE 2

Test(new_normal_random_vec, basic)
{
    MODULE* module = new_module_info_p(NBASE);
    double* pol_univ = malloc(NBASE * KBASE *sizeof(double));

    new_normal_random_vec(NBASE, pol_univ, KBASE, 2, 0.0, 0.001);

    for(int64_t p = 0 ; p < NBASE ; p++){
        printf("\n %e X^%ld \n", pol_univ[p], p);
    }

    delete_module_info_p(module);
    free(pol_univ);
}

Test(new_uniform_random_vec_dft, basic)
{
    MODULE* module = new_module_info_p(NBASE);
    PolyUnivDFT* res_dft = malloc(NBASE * KBASE *sizeof(double));

    new_uniform_random_vec_dft(module, res_dft, KBASE, 2);

    PolyUniv* res = malloc(NBASE * KBASE *sizeof(int64_t));

    vec_znx_dft_p(module, res_dft, KBASE, res, KBASE, NBASE);


    delete_module_info_p(module);
    free(res); free(res_dft);
}