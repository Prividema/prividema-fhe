#ifndef PVDA_BACKEND_NT
#define PVDA_BACKEND_NT

#include <stdint.h>

#include "backend.h"

/**
 * @file backend_nt.h
 *
 * @brief Number-theory functions in the Prividema-FHE backend
 *
 *
 */

/*
 * @brief Finds the primitive 2n-th root of ring Z_t
 *
 * @note t should be a prime that is 1 modulo 2*N
 *
 * @param backend The prividema backend
 */
uint64_t pvda_ring_2nth_root(const PvdaBackend* backend, uint64_t t);

#endif
