#ifndef HEINT_H
#define HEINT_H

#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"

int heint_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t t,
                 uint64_t* in);
int heint_decode(const PvdaBackend* backend, const GLWEParams* params, uint64_t* out, uint64_t slots, int64_t t,
                 PolyBiv* in);
#endif
