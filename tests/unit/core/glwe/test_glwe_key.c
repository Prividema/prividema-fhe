#include "core/glwe/glwe_key.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 4
#define KBASE 1
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE -12

//! GLWE KEY PART (begin)

/** 
 * @brief Ensures new_secret_key_values creates no NULL-pointer.
*/
Test(new_secret_key_values, basic){
    PolyUniv** values = new_secret_key_values(NBASE, KBASE);

    cr_assert(eq(int, values != NULL, 1));
    for(int64_t j = 0 ; j < KBASE ; j++)
        cr_assert(eq(int, values[j] != NULL, 1));

    delete_secret_key_values(values, KBASE);
}

/** 
 * @brief Ensure new_uniform_secret_key_values creates no NULL-pointer.
*/
Test(new_uniform_secret_key_values, basic){
    PolyUniv** values = new_uniform_secret_key_values(NBASE, KBASE, 3);

    cr_assert(eq(int, values != NULL, 1));
    for(int64_t j = 0 ; j < KBASE ; j++)
        cr_assert(eq(int, values[j] != NULL, 1));

    delete_secret_key_values(values, KBASE);
}

/**
 * @brief Tests wether transform_secret_key_values_dft_to_not_dft actually transforms the values of the secret key in DFT space, out of DFT space.
 */
Test(transform_secret_key_values_dft_to_not_dft, basic){
    PolyUnivDFT** values_dft = new_uniform_secret_key_values_dft(NBASE, KBASE, 3);
    PolyUniv** values = transform_secret_key_values_dft_to_not_dft(values_dft, NBASE, KBASE);

    for(int64_t j = 0 ; j < KBASE ; j++)
        for(int64_t p = 0 ; p < NBASE ; p++)
            cr_log_info("%" PRId64 " X^%" PRId64, values[j][p], p);

    delete_secret_key_values_dft(values_dft, KBASE);
    delete_secret_key_values(values, KBASE);
}

/**
 * @brief Ensures new_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key, values_not_null){
    PolyBiv** values = new_secret_key_values(KBASE, NBASE);
    GLWESecretKey* sk = new_glwe_secret_key(values, NBASE, KBASE);

    cr_assert(eq(int, sk != NULL, 1));
    cr_assert(eq(int, sk->values != NULL, 1));

    for(int64_t j = 0 ; j < KBASE ; j++)
        cr_assert(eq(int, sk->values[j] != NULL, 1));

    delete_glwe_secret_key(sk);
} 

/**
 * @brief Ensures new_uniform_glwe_secret_key returns a non-NULL pointer when values != NULL.
 */
Test(new_uniform_glwe_secret_key, values_not_null){
    GLWESecretKey* sk = new_uniform_glwe_secret_key(NBASE, KBASE, 2);

    cr_assert(eq(int, sk != NULL, 1));
    cr_assert(eq(int, sk->values != NULL, 1));

    for(int64_t j = 0 ; j < KBASE ; j++)
        cr_assert(eq(int, sk->values[j] != NULL, 1));

    delete_glwe_secret_key(sk);
} 

/**
 * @brief Tests wether transform_glwe_secret_key_dft_to_not_dft transforms the secret key in DFT space, out of DFT space.
 */
Test(transform_glwe_secret_key_dft_to_not_dft, basic){
    GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 3);
    GLWESecretKey* sk = transform_glwe_secret_key_dft_to_not_dft(sk_dft);

    for(int64_t j = 0 ; j < KBASE ; j++)
        for(int64_t p = 0 ; p < NBASE ; p++)
            cr_log_info("%" PRId64 " X^%" PRId64, sk->values[j][p], p);

    delete_glwe_secret_key_dft(sk_dft);
    delete_glwe_secret_key(sk);
}


//! GLWE KEY PART IN DFT SPACE (begin)


/** 
 * @brief Ensure new_secret_key_values_dft creates no NULL-pointer.
*/
Test(new_secret_key_values_dft, basic){
    PolyUnivDFT** values_dft = new_secret_key_values_dft(NBASE, KBASE);

    cr_assert(eq(int, values_dft != NULL, 1));
    for(int64_t j = 0 ; j < KBASE ; j++){
        cr_assert(eq(int, values_dft[j] != NULL, 1));
    }

    delete_secret_key_values_dft(values_dft, KBASE);
}

/** 
 * @brief Ensure new_uniform_secret_key_values_dft creates no NULL-pointer.
*/
Test(new_uniform_secret_key_values_dft, basic){
    PolyUnivDFT** values_dft = new_uniform_secret_key_values_dft(NBASE, KBASE, 3);

    cr_assert(eq(int, values_dft != NULL, 1));
    for(int64_t j = 0 ; j < KBASE ; j++){
        cr_assert(eq(int, values_dft[j] != NULL, 1));
    }

    delete_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Tests wether transform_secret_key_values_not_dft_to_dft actually transforms the values of the secret key out of DFT space, in DFT space.
 */
Test(transform_secret_key_values_not_dft_to_dft, basic){
    PolyUniv** values = new_uniform_secret_key_values(NBASE, KBASE, 3);
    PolyUnivDFT** values_dft = transform_secret_key_values_not_dft_to_dft(values, NBASE, KBASE);

    for(int64_t j = 0 ; j < KBASE ; j++){
        for(int64_t p = 0 ; p < NBASE ; p++){
            cr_log_info("%lf X^%" PRId64, values_dft[j][p], p);
        }
    }

    delete_secret_key_values(values, KBASE);
    delete_secret_key_values_dft(values_dft, KBASE);
}

/**
 * @brief Ensures new_glwe_secret_key_dft returns a non-NULL pointer when values != NULL.
 */
Test(new_glwe_secret_key_dft, values_not_null){
    PolyBivDFT** values = new_secret_key_values_dft(KBASE, NBASE);
    GLWESecretKeyDFT* sk_dft = new_glwe_secret_key_dft(values, NBASE, KBASE);

    cr_assert(eq(int, sk_dft != NULL, 1));
    cr_assert(eq(int, sk_dft->values != NULL, 1));

    for(int64_t j = 0 ; j < KBASE ; j++){
        cr_assert(eq(int, sk_dft->values[j] != NULL, 1));
    }

    delete_glwe_secret_key_dft(sk_dft);
} 

/**
 * @brief Tests whether new_uniform_glwe_secret_key_gen works as intended.
 */
Test(new_uniform_glwe_secret_key_dft, what_s_inside){
    GLWESecretKeyDFT* sk_dft = new_uniform_glwe_secret_key_dft(NBASE, KBASE, 2);

    cr_assert(eq(int, sk_dft != NULL, 1));
    cr_assert(eq(int, sk_dft->values != NULL, 1));

    for(int64_t j = 0 ; j < KBASE ; j++){
        cr_assert(eq(int, sk_dft->values[j] != NULL, 1));
    }
    
    delete_glwe_secret_key_dft(sk_dft);
}
     
/**
 * @brief Tests wether transform_glwe_secret_key_not_dft_to_dft transforms the secret key out of DFT space, in DFT space.
 */
Test(transform_glwe_secret_key_not_dft_to_dft, basic){
    GLWESecretKey* sk = new_uniform_glwe_secret_key(NBASE, KBASE, 3);
    GLWESecretKeyDFT* sk_dft = transform_glwe_secret_key_not_dft_to_dft(sk);

    for(int64_t j = 0 ; j < KBASE ; j++)
        for(int64_t p = 0 ; p < NBASE ; p++)
            cr_log_info("%lf X^%" PRId64, sk_dft->values[j][p], p);

    delete_glwe_secret_key(sk);
    delete_glwe_secret_key_dft(sk_dft);
}