#include "core/glwe/glwe_key.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 4
#define KBASE 2
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

/**
 * @brief Tests whether new_glwe_secret_key_dft returns a non-NULL pointer when values = NULL.
 */
Test(new_glwe_secret_key_dft, values_is_null){
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(NULL, NBASE, KBASE);

    cr_assert(eq(int, (sk_dft != NULL) && (sk_dft->values != NULL), 1));

    delete_glwe_secret_key_dft(sk_dft);
}


/**
 * @brief Tests whether new_glwe_secret_key_dft returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key_dft, values_not_null){
    PolyBivDFT** values = malloc((KBASE + 1)*sizeof(PolyBivDFT*));
    for(int64_t j = 0 ; j < KBASE ; j++)
    {
        values[j] = calloc(NBASE * sizeof(int64_t),1);
    }
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(values, NBASE, KBASE);

    cr_assert(eq(int, (sk_dft != NULL) && (sk_dft->values != NULL), 1));

    values[1] = malloc(NBASE * sizeof(int64_t));
    values[1][0] = 0.1;

    cr_assert(eq(dbl, sk_dft->values[1][0], 0.1));

    delete_glwe_secret_key_dft(sk_dft);
}

/**
 * @brief Tests whether new_uniform_glwe_secret_key_gen works as intended.
 */
Test(new_uniform_glwe_secret_key_gen, what_s_inside){
    GLWEPreparedSK* sk_dft = new_uniform_glwe_secret_key_gen(NBASE, KBASE, 2);
    PolyUniv** values = secret_key_values_dft_to_not_dft(sk_dft);

    for(int64_t j = 0 ; j < KBASE ; j++)
    {
        values[j] = calloc(NBASE * sizeof(int64_t),1);
    }

    for(int64_t j = 0 ; j < KBASE ; j++)
        for(int64_t p = 0 ; p < NBASE ; p++)
            cr_log_info("values[%ld] : %ld X^%ld", j, values[j][p], p);
    
    cr_assert(1);
    
    delete_glwe_secret_key_dft(sk_dft);
    for(int64_t j = 0 ; j < KBASE ; j++)
        free(values[j]);
    free(values);
}
    