#ifndef MONTARITH_H
#define MONTARITH_H

#include <stdint.h>
uint64_t montgomery_encode_32bit(uint64_t x, uint64_t q, uint64_t q_tild, uint64_t r2modq);
uint64_t montgomery_tild_32bit(uint64_t q);
uint64_t montgomery_r2modq_32bit(uint64_t q);
uint64_t montgomery_red_32bit(uint64_t m, uint64_t q, uint64_t q_tild);
//uint64_t montgomery_encode_32bit(uint64_t x, uint64_t q);
uint64_t montgomery_decode_32bit(uint64_t x, uint64_t q, uint64_t q_tild);
uint64_t montgomery_mult_32bit(uint64_t x_m, uint64_t y_m, uint64_t q, uint64_t q_tild);
uint64_t montgomery_pow_exp_32bit(uint64_t base_m, uint64_t exp, uint64_t q, uint64_t q_tild, uint64_t one_m);

#endif
