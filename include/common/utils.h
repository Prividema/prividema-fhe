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

#define CHECK_ALLOC_LABEL(ptr, m, label) \
	do                                   \
	{                                    \
		if (!(ptr))                      \
		{                                \
			log_perror(m);               \
			goto label;                  \
		}                                \
	} while (0)

/**
 * @brief Raises a runtime error
 *
 * @param m   The error message
 */
#define RAISE_ERROR(m)             \
	do                             \
	{                              \
		log_message(LOG_ERROR, m); \
		goto cleanup;              \
	} while (0)

/**
 * @brief Checks if an expression throws an error.
 *
 * This function checks if an expression (most of the time it's the return value of a function)
 * throws an error i.e its value is -1 (or a signed value < 0).
 *
 * @param expr The function to run.
 * @param m    The error message header.
 *
 * @note the `do ... while` isolates the `if` so we avoid errors if there's an `else` after a call to the macro.
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
 * @brief Checks if an expression throws an error, with programmable cleanup target
 *
 * This function checks if an expression (most of the time it's the return value of a function)
 * throws an error i.e its value is -1 (or a signed value < 0).
 *
 * @param expr The function to run.
 * @param m    The error message header.
 * @param label target label for the cleanup code
 *
 * @note the `do ... while` isolates the `if` so we avoid errors if there's an `else` after a call to the macro.
 */
#define CHECK_CALL_LABEL(expr, m, label) \
	do                                   \
	{                                    \
		if ((expr) < 0)                  \
		{                                \
			log_message(LOG_ERROR, m);   \
			goto label;                  \
		}                                \
	} while (0)

/**
 * @brief Computes the distance between two torus elements
 *
 * In other words, abs( (a mod 1) - (b mod 1))
 * Mostly used in tests
 *
 * @param a An element in the torus
 * @param b Another element in the torus
 *
 * @returns The (shortest) distance between the two elements in the torus
 */
double rnx_torus_distance(double a, double b);

/**
 * @brief Computes the distance between two torus elements
 *
 * In other words, abs( (a mod 1) - (b mod 1))
 * Mostly used in tests
 *
 * @param a An element in the torus
 * @param b Another element in the torus
 *
 * @returns The (shortest) distance between the two elements in the torus
 */
uint64_t tnx_torus_distance(uint64_t a, uint64_t b);

uint64_t u64_round_up_div(uint64_t num, uint64_t step);
int32_t i32_round_up_div(int32_t num, int32_t step);

/**
 * @brief Generic macro to compute an integer division with rounding up when the remainder is non-zero.
 *
 * @param num The dividend
 * @param step The divisior
 *
 */
#define INT_ROUND_UP_DIV(num, step) _Generic((num), uint64_t: u64_round_up_div, int32_t: i32_round_up_div)(num, step)

#endif  // UTILS_H
