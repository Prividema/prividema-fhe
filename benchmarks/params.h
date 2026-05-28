#ifndef PARAMS_H
#define PARAMS_H

#include <cmath>

#define NBASE      (1 << 14)
#define KBASE      1
#define KAPPABASE  19
#define NLIMBSBASE (15 * 2)
#define LBASE      NLIMBSBASE / (KBASE + 1)
#define SIGMABITS  4  //bits of sigma in the last limb. Should not affect performance
#define SIGMABASE  (ldexp(1.0, SIGMABITS - (LBASE) * KAPPABASE))
#define SKBITS     3  // bit size of secret key coefficients

#endif
