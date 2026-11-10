#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>
#include <string.h>
#include <criterion/redirect.h>

#include "errno.h"
#include "logger.h"

Test(test_one_log, test_log)
{
    printf("==== test_one_log ====\n");
    printf("Expecting 1 error message.\n");
    set_log_flag(LOG_ERROR);
    int a = log_message(LOG_DEBUG, "This is a  Debug message."); // Should be no print
    int b = log_message(LOG_INFO,  "This is an Info  message."); // Should be no print
    int c = log_message(LOG_WARN,  "This is a  Warn  message."); // Should be no print
    int d = log_message(LOG_ERROR, "This is an Error message.");

    cr_assert(a == 0);
    cr_assert(b == 0);
    cr_assert(c == 0);
    cr_assert(d < 0);
}

Test(test_all_logs, test_log)
{
    printf("==== test_all_logs ====\n");
    printf("Expecting 4 messages of different types.\n");
    set_log_flag(LOG_ERROR | LOG_WARN | LOG_INFO | LOG_DEBUG);
    int a = log_message(LOG_DEBUG, "This is a  Debug message.");
    int b = log_message(LOG_INFO,  "This is an Info  message.");
    int c = log_message(LOG_WARN,  "This is a  Warn  message.");
    int d = log_message(LOG_ERROR, "This is an Error message.");

    cr_assert(a == 0);
    cr_assert(b == 0);
    cr_assert(c == 0);
    cr_assert(d < 0);
    printf("\n");
}

Test(test_log_perror, test_log)
{
    printf("==== test_log_perror ====\n");
    printf("Expecting an error message.\n");
    
    errno = ENOMEM;
    int a = log_perror("Test ENOMEM");
    cr_assert(a < 0);
    printf("\n");
}

Test(test_log_is_null, test_log)
{
    printf("==== test_log_is_null ====\n");
    printf("Expecting a single error message.\n");
    int* p = NULL;
    errno = ENOMEM; // Simulates a malloc failure
    int np = 3;
    int a = log_is_null(p, "Test null pointer");
    int b = log_is_null(&np, "Test non-null pointer"); // Should be no print

    cr_assert(a < 0);
    cr_assert(b == 0);
    printf("\n");
}
Test(test_log_filtered, test_log)
{
    // A filtered message prints nothing, but an error is still reported as a failure
    set_log_flag(0);
    cr_assert(log_message(LOG_ERROR, "filtered error") < 0);
    cr_assert(log_message(LOG_WARN, "filtered warning") == 0);
    cr_assert(log_message(LOG_DEBUG, "filtered debug %d", 3) == 0);
    errno = EINVAL;
    cr_assert(log_perror("filtered perror") < 0);
    cr_assert(log_is_null(NULL, "filtered null") < 0);
}

Test(test_log_user_levels, test_log, .init = cr_redirect_stdout)
{
    set_log_flag(LOG_USR1 | LOG_USR4);
    cr_assert(log_message(LOG_USR1, "user %s", "one") == 0);
    cr_assert(log_message(LOG_USR2, "not printed") == 0);
    cr_assert(log_message(LOG_USR4, "user four") == 0);
    fflush(stdout);

    char out[256] = {0};
    FILE* f = cr_get_redirected_stdout();
    size_t len = fread(out, 1, sizeof(out) - 1, f);
    out[len] = '\0';
    cr_assert(strstr(out, "USR1: user one") != NULL, "output: %s", out);
    cr_assert(strstr(out, "USR4: user four") != NULL, "output: %s", out);
    cr_assert(strstr(out, "not printed") == NULL, "output: %s", out);
}

Test(test_log_format, test_log, .init = cr_redirect_stdout)
{
    set_log_flag(LOG_INFO);
    cr_assert(log_message(LOG_INFO, "value=%d text=%s", 42, "ok") == 0);
    fflush(stdout);

    char out[256] = {0};
    FILE* f = cr_get_redirected_stdout();
    size_t len = fread(out, 1, sizeof(out) - 1, f);
    out[len] = '\0';
    // "[YYYY-MM-DD HH:MM:SS] INFO : value=42 text=ok\n"
    cr_assert(out[0] == '[' && out[20] == ']', "output: %s", out);
    cr_assert(strstr(out, "] INFO : value=42 text=ok\n") != NULL, "output: %s", out);
}
