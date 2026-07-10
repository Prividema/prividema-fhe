#ifndef CKKS_H
#define CKKS_H

#include <complex.h>
#include <stdint.h>

#include "backend.h"

void encode_slow_internal(double* out, uint64_t n, double complex* in);
void complex_bitrev(double complex* inout, uint64_t n);
int encode_internal(const PvdaBackend* backend, double* out, uint64_t slots, double complex* in, int inplace);
void decode_slow_internal(double complex* out, uint64_t n, double* in);

#endif
