#ifndef UTILS_H
#define UTILS_H

#include "logger.h"
#include "stdint.h"

#endif  // UTILS_H
#define CHECK_ALLOC(ptr, m) \
	do                        \
	{                         \
		if (!(ptr))           \
		{                     \
			log_perror(m);  \
			goto cleanup;     \
		}                     \
	} while (0)

#define CHECK_CALL(expr, m) \
	do                        \
	{                         \
		if ((expr) < 0)       \
		{                     \
			log_perror(m);  \
			goto cleanup;     \
		}                     \
	} while (0)



/**
 * @brief Returns the probality of \f$ X \leq \mathbb{E}(X) \f$, where X is the random variable that counts the number of error greater or equal than 3*sigma.
 * 
 * @note The errors are drawn normally, the standard deviation is sigma.
 * @param N The number of draws.
 * @param sigma The standard deviation of the normal error's distribution.
 * @return double 
 */
double binomial_tail(uint64_t N, double p);