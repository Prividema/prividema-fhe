#ifndef UTILS_H
#define UTILS_H

#endif  // UTILS_H
#define CHECK_ALLOC(ptr, msg) \
	do                        \
	{                         \
		if (!(ptr))           \
		{                     \
			log_perror(msg);  \
			goto cleanup;     \
		}                     \
	} while (0)

#define CHECK_CALL(expr, msg) \
	do                        \
	{                         \
		if ((expr) < 0)       \
		{                     \
			log_perror(msg);  \
			goto cleanup;     \
		}                     \
	} while (0)
