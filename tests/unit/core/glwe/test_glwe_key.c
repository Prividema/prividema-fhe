#include "glwe_key.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#define NBASE 1024
#define KBASE 8
#define KAPPABASE 4
#define NLIMBSBASE 45
#define LBASE NLIMBSBASE/(KBASE+1)
#define SIGMABASE 1e-7

Test(new_glwe_secret_key_dft, values_is_null){
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(NULL, NBASE, KBASE);

    cr_assert(eq(int, (sk_dft != NULL) && (sk_dft->values != NULL), 1));

    delete_glwe_secret_key_dft(sk_dft);
}

// TODO
Test(new_glwe_secret_key_dft, values_not_null){
    PolyBivDFT** values = malloc((KBASE + 1)*sizeof(PolyBivDFT*));
    GLWEPreparedSK* sk_dft = new_glwe_secret_key_dft(values, NBASE, KBASE);
}

/**
 * @brief Creates a GLWE Secret key in DFT space
 * 
 * @param values The values of the secret key in DFT space.
 * @param N The polynomials' maximum degree in X.
 * @param k The number of ZnX polynomial in the secret key.
 * @return GLWEPreparedSK* 
 */
GLWEPreparedSK* new_glwe_secret_key_dft(PolyUnivDFT** values, int64_t N, int64_t k);

/**
 * @brief Delete the secret key that is in DFT space.
 * 
 * @param sk_dft The secret key in DFT space.
 */
void delete_glwe_secret_key_dft(GLWEPreparedSK* sk_dft);