#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <errno.h>
#include <stdlib.h>

#include "logger.h"
#include "utils.h"

// Tests of the error-handling macros of utils.h: each one must jump to its label exactly when the condition fails,
// and fall through otherwise.

static int check_alloc(void* ptr, int* reached_after)
{
	*reached_after = 0;
	CHECK_ALLOC(ptr, "test allocation");
	*reached_after = 1;
	return 0;
cleanup:
	return -1;
}

static int check_alloc_label(void* ptr)
{
	CHECK_ALLOC_LABEL(ptr, "test allocation", failed);
	return 0;
failed:
	return -1;
}

static int check_call(int value, int* reached_after)
{
	*reached_after = 0;
	CHECK_CALL(value, "test call");
	*reached_after = 1;
	return 0;
cleanup:
	return -1;
}

static int check_call_label(int value)
{
	CHECK_CALL_LABEL(value, "test call", failed);
	return 0;
failed:
	return -1;
}

static int raise_error(int raise)
{
	if (raise) RAISE_ERROR("raised");
	return 0;
cleanup:
	return -1;
}

static int raise_error_label(int raise)
{
	if (raise) RAISE_ERROR_LABEL("raised", failed);
	return 0;
failed:
	return -1;
}

/** @brief The macros are single statements: they can be used before an else without braces. */
static int macro_before_else(int value)
{
	int r = 0;
	if (value > 0)
		CHECK_CALL(value - 10, "negative in branch");
	else
		r = 2;
	return r;
cleanup:
	return -1;
}

Test(utils_macros, check_alloc)
{
	set_log_flag(0);  // keep the expected errors quiet
	int x = 0, reached;
	errno = ENOMEM;
	cr_assert(eq(int, check_alloc(NULL, &reached), -1));
	cr_assert(eq(int, reached, 0));
	cr_assert(eq(int, check_alloc(&x, &reached), 0));
	cr_assert(eq(int, reached, 1));
	cr_assert(eq(int, check_alloc_label(NULL), -1));
	cr_assert(eq(int, check_alloc_label(&x), 0));
}

Test(utils_macros, check_call)
{
	set_log_flag(0);
	int reached;
	cr_assert(eq(int, check_call(-1, &reached), -1));
	cr_assert(eq(int, reached, 0));
	cr_assert(eq(int, check_call(-42, &reached), -1));  // any negative value is an error
	cr_assert(eq(int, check_call(0, &reached), 0));
	cr_assert(eq(int, reached, 1));
	cr_assert(eq(int, check_call(5, &reached), 0));  // positive values are successes
	cr_assert(eq(int, check_call_label(-1), -1));
	cr_assert(eq(int, check_call_label(0), 0));
}

Test(utils_macros, raise_error)
{
	set_log_flag(0);
	cr_assert(eq(int, raise_error(1), -1));
	cr_assert(eq(int, raise_error(0), 0));
	cr_assert(eq(int, raise_error_label(1), -1));
	cr_assert(eq(int, raise_error_label(0), 0));
}

Test(utils_macros, single_statement)
{
	set_log_flag(0);
	cr_assert(eq(int, macro_before_else(20), 0));
	cr_assert(eq(int, macro_before_else(5), -1));
	cr_assert(eq(int, macro_before_else(-1), 2));
}
