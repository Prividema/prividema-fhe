#include "common/utils.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>

// Test gsl_init
Test(common, test_gsl_rng_init) 
{
    gsl_rng *r;
    gsl_init(r);
    cr_assert(ne(r, NULL), "gsl_init produced a NULL pointer. It should be non NULL");
    gsl_rng_free(r);
}

// Test gsl_refresh_seed
Test(common, test_gsl_refresh_seed) 
{
    gsl_rng *r;
    gsl_init(r);
    gsl_refresh_seed(r);
    int rng1 = gsl_rng_uniform_int(r, UINT32_MAX);
    gsl_refresh_seed(r);
    int rng2 = gsl_rng_uniform_int(r, UINT32_MAX);
    cr_assert_neq(rng1, rng2, "Two different seeds generated the same first random number. You are either not changing the seed or extremely unlucky !");
    gsl_rng_free(r);
}