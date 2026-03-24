#ifndef UTILS_H
#define UTILS_H

#include <stdint.h>

#include "logger.h"

/**
 * @brief Checks if a pointer is allocated.
 *
 * @param ptr The pointer to test.
 * @param m	  The perror message header.
 *
 * @note the `do ... while` isolates the `if`
 * 		 so we avoid errors if there's an `else` after a call to the macro.
 */
#define CHECK_ALLOC(ptr, m) \
	do                      \
	{                       \
		if (!(ptr))         \
		{                   \
			log_perror(m);  \
			goto cleanup;   \
		}                   \
	} while (0)

/**
 * @brief Checks if an expression throws an error.
 *
 * This function checks if an expression (most of the time it's the return value of a function)
 * throws an error i.e its value is -1 (or a signed value < 0).
 *
 * @param expr The function to run.
 * @param m	   The error message header.
 *
 * @note the `do ... while` isolates the `if`
 * 		 so we avoid errors if there's an `else` after a call to the macro.
 */
#define CHECK_CALL(expr, m)            \
	do                                 \
	{                                  \
		if ((expr) < 0)                \
		{                              \
			log_message(LOG_ERROR, m); \
			goto cleanup;              \
		}                              \
	} while (0)

/**
 * @brief Computes the distance between two torus elements
 *
 * In other words, abs( (a mod 1) - (b mod 1))
 * Mostly used in tests
 */
double torus_distance(double a, double b);

double binomial_tail(uint64_t nn, double p, int k);
#endif  // UTILS_H
