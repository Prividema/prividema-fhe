#ifndef CKKS_H
#define CKKS_H

#include <complex.h>
#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"

void encode_slow_internal(double* out, uint64_t n, double complex* in);
void decode_slow_internal(double complex* out, uint64_t n, double* in);
int encode_internal(const PvdaBackend* backend, double* out, uint64_t slots, double complex* in, int inplace);
int decode_internal(const PvdaBackend* backend, double complex* out, uint64_t slots, double* in);

int ckks_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t scale_bits,
                double complex* in);
int ckks_decode(const PvdaBackend* backend, const GLWEParams* params, double complex* out, uint64_t slots,
                int64_t scale_bits, PolyBiv* in);
#endif
