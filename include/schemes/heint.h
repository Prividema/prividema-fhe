#ifndef HEINT_H
#define HEINT_H

#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"

uint64_t mod_inv(int64_t x, int64_t mod);
uint64_t montgomery_red_32bit(uint64_t m, uint64_t q, uint64_t q_tild);
uint64_t montgomery_encode_32bit(uint64_t x, uint64_t q);
uint64_t montgomery_decode_32bit(uint64_t x, uint64_t q, uint64_t q_tild);
uint64_t montgomery_mult_32bit(uint64_t x_m, uint64_t y_m, uint64_t q, uint64_t q_tild);
uint64_t montgomery_pow_exp_32bit(uint64_t base_m, uint64_t exp, uint64_t q, uint64_t q_tild, uint64_t one_m);

int heint_encode(const PvdaBackend* backend, const GLWEParams* params, PolyBiv* out, uint64_t slots, int64_t t,
                 uint64_t* in);
int heint_decode(const PvdaBackend* backend, const GLWEParams* params, uint64_t* out, uint64_t slots, int64_t t,
                 PolyBiv* in);

int internal_slow_intt_heint(const PvdaBackend* backend, uint64_t nn, uint64_t* root_table_m, uint64_t* out_int,
                             uint64_t* in, uint64_t t);
int internal_slow_ntt_heint(const PvdaBackend* backend, uint64_t nn, uint64_t* root_table_m, uint64_t* out_int,
                            uint64_t* in, uint64_t t);
#endif
