#include "core/glwe/glwe_ciphertext.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

/**
 * @brief Test new_glwe
 */
Test(new_glwe, basic){
    GLWECtParams* params = new_glwe_ct_params(NBASE, KBASE, KAPPABASE, NLIMBSBASE, SIGMABASE);
    GLWECiphertext* ct = new_glwe(params);
    
    cr_assert(eq(int, (ct != NULL)&&(ct->vec != NULL), 1));

    delete_glwe(ct);
    delete_glwe_ct_params(params);
}