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
typedef struct pvda_flat_poly_vec_t
{
    uint64_t nn;    ///< Ring dimension (length of each polynomial)
    uint64_t size;  ///< Number of Znx polynomials
    int64_t stride; ///< Distance between consecutive polynomials, in elements (nn when contiguous)
    int64_t* ptr;   ///< Pointer to the first coefficient (not owned)
} FlatPolyVec;

static inline FlatPolyVec flat_poly_vec_view(uint64_t nn, uint64_t size, int64_t stride, int64_t* ptr)
{
    FlatPolyVec v = {nn,size, stride, ptr};
    return v;
}

// Opaque types

typedef struct pvda_be_poly_univ_dft_t      BE_PolyUnivZnxDFT;  ///< Backend univariate polynomial in DFT form.
typedef struct pvda_be_poly_biv_dft_t       BE_PolyBivDFT;      ///< Backend bivariate polynomials in DFT form.

typedef struct pvda_be_poly_univ_svp_prep_t BE_PolyUnivSvpPrep; ///< Backend univariate polynomial(scalar) prepared for scalar-vector product.
// In spqlios, the vector in svp is not prepared so there's no need for a backend structure for it, just a FlatPolyVec is sufficient.

typedef struct pvda_be_poly_vec_vmp_prep_t  BE_PolyBivVmpPrep;   ///< Backend vector of polynomials prepared for vector-matrix product
typedef struct pvda_be_poly_mat_vmp_prep_t  BE_PolyMatVmpPrep;   ///< Backend matrix of polynomials prepared for vector-matrix product


// Zero-cost: same memory, only the type changes.
// A univariate DFT polynomial has the same layout as a bivariate DFT polynomial with 1 limb
static inline BE_PolyBivDFT* univ_znx_dft_as_biv(BE_PolyUnivZnxDFT* p) { return (BE_PolyBivDFT*)p; }
static inline const BE_PolyBivDFT* univ_znx_dft_as_biv_const(const BE_PolyUnivZnxDFT* p) { return (const BE_PolyBivDFT*)p; }

#endif // PVDA_BACKEND_TYPES_H