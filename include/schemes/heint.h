#ifndef HEINT_H
#define HEINT_H

#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"

int heint_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t t,
                 uint64_t* in);
int heint_decode(const PvdaBackend* backend, const GLWEParams* params, uint64_t* out, uint64_t slots, int64_t t,
                 PolyBiv* in);

int internal_slow_intt_heint(const PvdaBackend* backend, uint64_t nn, uint64_t* root_table_m, uint64_t* out_int,
                             uint64_t* in, uint64_t t);
int internal_slow_ntt_heint(const PvdaBackend* backend, uint64_t nn, uint64_t* root_table_m, uint64_t* out_int,
                            uint64_t* in, uint64_t t);
#endif
