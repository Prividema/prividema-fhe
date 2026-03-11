#ifndef UTILS_H
#define UTILS_H

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
