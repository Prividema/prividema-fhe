#include <stdio.h>
#include <criterion/criterion.h>
#include "glwe.h"

// Test the add function
Test(math_utils, test_add) {
    cr_assert_eq(add(2, 3), 5, "2 + 3 should equal 5");
    cr_assert_eq(add(-1, 1), 0, "-1 + 1 should equal 0");
}

// Test the multiply function
Test(math_utils, test_multiply) {
    cr_assert_eq(multiply(4, 5), 20, "4 * 5 should equal 20");
    cr_assert_eq(multiply(-2, 3), -6, "-2 * 3 should equal -6");
}
