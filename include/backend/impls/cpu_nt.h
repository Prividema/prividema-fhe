#ifndef CPU_NT_H
#define CPU_NT_H
#include <stdint.h>

uint64_t pow_exp(uint64_t base, uint64_t exp, uint64_t m);
void pvda_fill_cpu_nt(struct pvda_virtual_table* vt);
#endif
