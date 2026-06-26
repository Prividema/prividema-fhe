#ifndef CKKS_H
#define CKKS_H

#include <complex.h>
#include <stdint.h>

void encode_slow_internal(double* out, uint64_t n, double complex* in);
void decode_slow_internal(double complex* out, uint64_t n, double* in);

#endif
