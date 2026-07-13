#ifndef HEFP_H
#define HEFP_H

#ifndef _cplusplus
#include <complex.h>
#endif
#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"

void hefp_encode_slow_internal(double* out, uint64_t n, double _Complex* in);
void hefp_decode_slow_internal(double _Complex* out, uint64_t n, double* in);
int hefp_encode_internal(const PvdaBackend* backend, double* out, uint64_t slots, double _Complex* in, int inplace);
int hefp_decode_internal(const PvdaBackend* backend, double _Complex* out, uint64_t slots, double* in);

int hefp_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t scale_bits,
                double _Complex* in);
int hefp_decode(const PvdaBackend* backend, const GLWEParams* params, double _Complex* out, uint64_t slots,
                int64_t scale_bits, PolyBiv* in);
#endif
