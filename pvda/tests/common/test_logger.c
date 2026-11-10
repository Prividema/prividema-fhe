#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <stdio.h>

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