#ifndef CPU_NT_H
#define CPU_NT_H
#include <stdint.h>

/**
 * @file cpu_nt.h
 *
 * CPU-based implementation of Number-Theory operations
 *
 */

/**
 * @brief Computes base^exp mod m
 *
 * @warning Non-constant-time on exp
 *
 * Consider using the montgomery alternative montgomery_pow_exp_32bit
 * for better performance
 *
 * @param base The base
 * @param exp The exmponent
 * @param m the modulus
 *
 * @return base^exp mod m
 */
uint64_t pow_exp(uint64_t base, uint64_t exp, uint64_t m);

/**
 * @brief Popluate virtual table enctries for CPU NT implementation
 *
 * @param vt The prividema virtual table to populate
 */
void pvda_fill_cpu_nt(struct pvda_virtual_table* vt);

#endif
