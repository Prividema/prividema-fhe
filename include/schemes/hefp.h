#ifndef HEFP_H
#define HEFP_H

#ifndef _cplusplus
#include <complex.h>
#endif
#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"

int hefp_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t scale_bits,
                double _Complex* in);
int hefp_decode(const PvdaBackend* backend, const GLWEParams* params, double _Complex* out, uint64_t slots,
                int64_t scale_bits, PolyBiv* in);
#endif
