#include "core/glwe/glwe.h"
#include "core/glwe/glwe_ct_params.h"
#include "vec_znx_arithmetic_private.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <time.h>
#include <math.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 32
#define NLIMBSBASE 180
#define LBASE NLIMBSBASE/(KBASE+1)


Test(vec_znx_normalize_base2k_p, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE);
    MODULE* module = new_module_info(params->N, FFT64);
    
    PolyBiv* a = malloc(LBASE * NBASE * sizeof(int64_t));
    PolyBiv* res = malloc(LBASE * NBASE * sizeof(int64_t));
    vec_znx_normalize_base2k_p(module, KAPPABASE, res, LBASE, NBASE, a, LBASE, NBASE);

    free(a); free(res);
    delete_glwe_ct_params(params);
    delete_module_info(module);
}