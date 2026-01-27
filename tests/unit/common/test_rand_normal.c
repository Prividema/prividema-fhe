#include "common/distributions.h"

#include <criterion/criterion.h>
#include <criterion/new/assert.h>

#include <stdio.h> 

#define NBASE 4
#define KBASE 2

Test(rand_normal, basic)
{
    double* a = malloc(sizeof(double));

    rand_normal(a, 0.0, 0.001);

    printf("%lf\n", *a);

    free(a);
}