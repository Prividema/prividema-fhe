#include "backend.h"

#include <stdlib.h>

#include "backend_private.h"
#include "impls/cpu_nt.h"
#include "impls/rng_cpu.h"
#include "impls/spqlios_alias.h"
#include "pvda_ffts.h"
#include "rng.h"

#ifdef OPENSSL_BACKEND
#include "impls/rng_openssl.h"
#endif

PvdaBackend* pvda_new_spqlios_backend(int nn)
{
	PvdaBackend* res    = malloc(sizeof(PvdaBackend));
	res->spqlios_module = spqlios_new_module_info(nn);
	pvda_fill_spqlios(&res->vt);
#ifdef OPENSSL_BACKEND
	pvda_fill_openssl_rng(&res->vt);
#else
	pvda_fill_ref_rng(&res->vt);
#endif
	res->pvda_fft_data = new_fft_data(nn);
	pvda_fill_cpu_nt(&res->vt);
	return res;
}

void pvda_delete_backend(PvdaBackend* module)
{
	spqlios_delete_module_info(module);
	delete_fft_data(module->pvda_fft_data);
	free(module);
}

// =============================================
// |                                           |
// |              BE_PolyUnivZnxDFT            |
// |                                           |
// =============================================

BE_PolyUnivZnxDFT* pvda_new_univ_znx_dft(const PvdaBackend* module)
{
	return module->vt.pvda_new_univ_znx_dft(module);
}

void pvda_delete_univ_znx_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol)
{
	module->vt.pvda_delete_univ_znx_dft(module, pol);
}

int pvda_flat_univ_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res, const FlatPolyVec* a)
{
	return module->vt.pvda_flat_univ_to_dft(module, res, a);
}

int pvda_flat_univ_idft(const PvdaBackend* module, FlatPolyVec* res, const BE_PolyUnivZnxDFT* pol)
{
	return module->vt.pvda_flat_univ_idft(module, res, pol);
}

// =============================================
// |                                           |
// |             BE_PolyUnivSvpPrep            |
// |                                           |
// =============================================

BE_PolyUnivSvpPrep* pvda_new_univ_svp_prep(const PvdaBackend* module)
{
	return module->vt.pvda_new_univ_svp_prep(module);
}

void pvda_delete_univ_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* pol)
{
	module->vt.pvda_delete_univ_svp_prep(module, pol);
}

int pvda_flat_univ_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* res, const FlatPolyVec *a)
{
	return module->vt.pvda_flat_univ_to_svp_prep(module, res, a);
}

// =============================================
// |                                           |
// |                BE_PolyBivDFT              |
// |                                           |
// =============================================

BE_PolyBivDFT* pvda_new_biv_dft(const PvdaBackend* module, uint64_t size)
{
	return module->vt.pvda_new_biv_dft(module, size);
}
void pvda_delete_biv_dft(const PvdaBackend* module, BE_PolyBivDFT* res)
{
	module->vt.pvda_delete_biv_dft(module, res);
}

int pvda_flat_biv_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const FlatPolyVec* a)
{
	return module->vt.pvda_flat_biv_to_dft(module, res, res_size, a);
}

int pvda_flat_biv_idft(const PvdaBackend* module, FlatPolyVec* res, const BE_PolyBivDFT* pol, uint64_t pol_size)
{
	return module->vt.pvda_flat_biv_idft(module, res, pol, pol_size);
}

// One vtable slot (const signature) serves both versions: the non-const one only casts back
// the pointer that the caller passed as non-const.
BE_PolyBivDFT* pvda_biv_dft_at(const PvdaBackend* module, BE_PolyBivDFT* vec, uint64_t offset)
{
	return (BE_PolyBivDFT*)module->vt.pvda_biv_dft_at(module, vec, offset);
}

const BE_PolyBivDFT* pvda_biv_dft_at_const(const PvdaBackend* module, const BE_PolyBivDFT* vec, uint64_t offset)
{
	return module->vt.pvda_biv_dft_at(module, vec, offset);
}

// =============================================
// |                                           |
// |                    SVP                    |
// |                                           |
// =============================================

int pvda_svp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                       const BE_PolyUnivSvpPrep* prepared_pol, const FlatPolyVec* a)
{
	return module->vt.pvda_svp_apply_dft(module, res, res_size, prepared_pol, a);
}

int pvda_svp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size,
                              const BE_PolyUnivSvpPrep* ppol, const BE_PolyBivDFT* a, uint64_t a_size)
{
	return module->vt.pvda_svp_apply_dft_to_dft(module, res, res_size, ppol, a, a_size);
}

// =============================================
// |                                           |
// |             BE_PolyBivVmpPrep             |
// |                                           |
// =============================================

BE_PolyBivVmpPrep* pvda_new_biv_vmp_prep(const PvdaBackend* module, uint64_t n_rows)
{
	return module->vt.pvda_new_biv_vmp_prep(module, n_rows);
}
void pvda_delete_biv_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* biv_prep)
{
	module->vt.pvda_delete_biv_vmp_prep(module, biv_prep);
}

int pvda_flat_biv_to_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* res, uint64_t res_size, const FlatPolyVec* a)
{
	return module->vt.pvda_flat_biv_to_vmp_prep(module, res, res_size, a);
}

// =============================================
// |                                           |
// |              BE_PolyMatVmpPrep            |
// |                                           |
// =============================================

BE_PolyMatVmpPrep* pvda_new_vmp_pmat(const PvdaBackend* module, uint64_t nrows, uint64_t ncols)
{
	return module->vt.pvda_new_vmp_pmat(module, nrows, ncols);
}

int pvda_vmp_prepare_contiguous(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat, const int64_t* mat, uint64_t nrows, uint64_t ncols)
{
	return module->vt.pvda_vmp_prepare_contiguous(module, pmat, mat, nrows, ncols);
}

void pvda_delete_vmp_pmat(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat)
{
	module->vt.pvda_delete_vmp_pmat(module, pmat);
}

// =============================================
// |                                           |
// |                    VMP                    |
// |                                           |
// =============================================

int pvda_vmp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const FlatPolyVec* a,
                       const BE_PolyMatVmpPrep* pmat, uint64_t nrows, uint64_t ncols)
{
	return module->vt.pvda_vmp_apply_dft(module, res, res_size, a, pmat, nrows, ncols);
}

int pvda_vmp_apply_dft_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                              const BE_PolyBivDFT* a_dft, uint64_t a_size, const BE_PolyMatVmpPrep* pmat,
                              const uint64_t nrows, const uint64_t ncols)
{
	return module->vt.pvda_vmp_apply_dft_to_dft(module, res, res_size, a_dft, a_size, pmat, nrows, ncols);
}

int pvda_vmp_apply_prepared_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, const uint64_t res_size,
                                   const BE_PolyBivVmpPrep* a_prep, uint64_t a_size, const BE_PolyMatVmpPrep* pmat,
                                   const uint64_t nrows, const uint64_t ncols)
{
	return module->vt.pvda_vmp_apply_prepared_to_dft(module, res, res_size, a_prep, a_size, pmat, nrows, ncols);
}

// =============================================
// |                                           |
// |                 Arithmetic                |
// |                                           |
// =============================================

int pvda_vec_znx_normalize_base2k(const PvdaBackend* module, uint64_t log2_base2k, FlatPolyVec* res, const FlatPolyVec* a)
{
	return module->vt.pvda_vec_znx_normalize_base2k(module, log2_base2k, res, a);
}

int pvda_znx_small_product(const PvdaBackend* module, int64_t* res, const int64_t* a, const int64_t* b)
{
	return module->vt.pvda_znx_small_product(module, res, a, b);
}

int pvda_vec_znx_negate(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a)
{
	return module->vt.pvda_vec_znx_negate(module, res, a);
}

int pvda_vec_znx_add(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a, const FlatPolyVec* b)
{
	return module->vt.pvda_vec_znx_add(module, res, a, b);
}

int pvda_vec_znx_sub(const PvdaBackend* module, FlatPolyVec* res, const FlatPolyVec* a, const FlatPolyVec* b)
{
	return module->vt.pvda_vec_znx_sub(module, res, a, b);
}

int pvda_znx_automorphism(const PvdaBackend* module, const int64_t p, int64_t* res, const int64_t* a)
{
	return module->vt.pvda_znx_automorphism(module, p, res, a);
}

int pvda_vec_znx_automorphism(const PvdaBackend* module, const int64_t p, FlatPolyVec* res, const FlatPolyVec* a)
{
	return module->vt.pvda_vec_znx_automorphism(module, p, res, a);
}

int pvda_vec_znx_rotate(const PvdaBackend* module, const int64_t p, FlatPolyVec* res, const FlatPolyVec* a)
{
	return module->vt.pvda_vec_znx_rotate(module, p, res, a);
}

uint64_t pvda_module_extract_nn(const PvdaBackend* module)
{
	return module->vt.pvda_module_extract_nn(module);
}

int pvda_rand_uniform_pow2(const PvdaBackend* module, int64_t* res, uint64_t nb_bits)
{
	return module->vt.pvda_rand_uniform_pow2(module, res, nb_bits);
}

int pvda_rand_uniform_pow2_vec(const PvdaBackend* module, int64_t* res, uint64_t n, uint64_t nb_bits)
{
	return module->vt.pvda_rand_uniform_pow2_vec(module, res, n, nb_bits);
}

int pvda_rand_uniform_binary_vec(const PvdaBackend* module, uint64_t* res, uint64_t n)
{
	return module->vt.pvda_rand_uniform_binary_vec(module, res, n);
}

int pvda_rand_uniform(const PvdaBackend* module, int64_t* res, int64_t limit_down, int64_t limit_up)
{
	return module->vt.pvda_rand_uniform(module, res, limit_down, limit_up);
}

uint64_t pvda_ring_2nth_root(const PvdaBackend* module, uint64_t t)
{
	return module->vt.pvda_ring_2nth_root(module, t);
}

int pvda_vec_znx_dft_add(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_n_limbs, BE_PolyBivDFT* a,
                         uint64_t a_n_limbs, BE_PolyBivDFT* b, uint64_t b_n_limbs)
{
	return module->vt.pvda_vec_znx_dft_add(module, res, res_n_limbs, a, a_n_limbs, b, b_n_limbs);
}