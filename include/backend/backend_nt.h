#ifndef PVDA_BACKEND_NT
#define PVDA_BACKEND_NT

#include <stdint.h>

#include "backend.h"

/*
 * @brief Finds the primitive 2n-th root of ring Z_t
 *
 *
 */
uint64_t pvda_ring_2nth_root(const PvdaBackend* backend, uint64_t t);

#endif
