#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Generates a random seed.
 * 
 * @param seed A pointer to the seed that has to be filled.
 *
 * @retval - `-1` if an error occurs. In this case, the error comes from the syscalls, and perror is called.
 * @retval - `0` otherwise.
 */
int rand_seed(unsigned long *seed)
{
    FILE *f = fopen("/dev/urandom", "rb");

    if (!f) {
        perror("fopen");
        return -1;
    }

    if (fread(seed, sizeof(seed), 1, f) != 1) {
        perror("fread");
        fclose(f);
        return -1;
    }
    fclose(f);
    return 0;
}