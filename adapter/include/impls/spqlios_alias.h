#ifndef SPQLIOS_ALIAS_H
#define SPQLIOS_ALIAS_H

#include <stdint.h>

#include "backend.h"
#include "backend_private.h"
#include "pvda_types.h"

/**
 *
 * @file spqlios_alias.h
 *
 * Adapter between prividema and the spqlios-arithmetic backend: creation and deletion of a spqlios module,
 * conversion of prividema data into spqlios input data, and registration of the spqlios functions in the
 * prividema virtual table.
 *
 */

// =============================================
// |                                           |
// |        spqlios coefficient-domain input   |
// |                                           |
// =============================================

/**
 * @brief Flat array of `size` Znx polynomials in host memory: the coefficient-domain exchange format of spqlios.
 *
 * Polynomial i occupies ptr[i * stride .. i * stride + nn - 1]. This matches the (ptr, size, stride) triplets that
 * the spqlios vec_znx functions take. It is specific to spqlios: other backends have their own input formats.
 */
typedef struct pvda_flat_poly_vec_t
{
	uint64_t nn;     ///< Ring dimension (length of each polynomial)
	uint64_t size;   ///< Number of Znx polynomials
	int64_t stride;  ///< Distance between consecutive polynomials, in elements (nn when contiguous)
	int64_t* ptr;    ///< Pointer to the first coefficient (not owned)
} FlatPolyVec;

/** @brief spqlios view of the limbs of a prividema bivariate polynomial. Does not copy. */
static inline FlatPolyVec spqlios_flat_from_biv(const PolyBiv* a)
{
	FlatPolyVec v = {a->nn, a->nb_limbs, a->stride, a->ptr};
	return v;
}

/** @brief spqlios view of a prividema univariate polynomial of nn coefficients (one polynomial). Does not copy. */
static inline FlatPolyVec spqlios_flat_from_univ(uint64_t nn, const PolyUnivZnx* a)
{
	FlatPolyVec v = {nn, 1, (int64_t)nn, (int64_t*)a};
	return v;
}

// =============================================
// |                                           |
// |                 Module                    |
// |                                           |
// =============================================

SPQLIOS_MODULE* spqlios_new_module_info(uint64_t nn);

void spqlios_delete_module_info(PvdaBackend* module);

void pvda_fill_spqlios(struct pvda_virtual_table* vt);
#endif
