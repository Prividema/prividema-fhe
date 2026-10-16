#ifndef PVDA_BACKEND_TYPES_H
#define PVDA_BACKEND_TYPES_H

#include <stdint.h>

/**
 * @file backend_types.h
 * @brief Types exchanged between the library and the backend
 */

/**
 * @brief Flat array of `size` Znx polynomials in host memory : the coefficient-domain exchange format with the backend.
 * 
 * Polynomial i occupies ptr[i * stride .. i * stride + nn - 1].
 */
typedef struct pvda_znx_vec_t
{
    uint64_t nn;    ///< Ring dimension (length of each polynomial)
    uint64_t size;  ///< Number of Znx polynomials
    int64_t stride; ///< Distance between consecutive polynomials, in elements (nn when contiguous)
    int64_t* ptr;   ///< Pointer to the first coefficient (not owned)
} ZnxVec;

static inline ZnxVec znx_vec_view(uint64_t nn, uint64_t size, int64_t stride, int64_t* ptr)
{
    ZnxVec v = {nn,size, stride, ptr};
    return v;
}

#endif // PVDA_BACKEND_TYPES_H