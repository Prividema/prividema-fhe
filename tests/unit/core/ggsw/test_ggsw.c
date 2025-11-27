#include "core/ggsw/ggsw.h"
#include <criterion/criterion.h>
#include <criterion/new/assert.h>

// Test the add function
Test(ggsw, test_add) {
    int a[1] = {1};
    int b[1] = {1};
    int c[1] = {2};
    cr_assert(eq(int[1],add(a,1,b,1), c), "add(2, 3) should be 5 but was %d", add(a,1,b,1));
}

// Test the multiply function
Test(math_utils, test_multiply) {
    cr_assert_eq(multiply(4, 5), 20, "4 * 5 should equal 20");
    cr_assert_eq(multiply(-2, 3), -6, "-2 * 3 should equal -6");
}
