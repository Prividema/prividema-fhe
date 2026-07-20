#include "backend.h"

#include <stdlib.h>

#include "backend_private.h"
#include "impls/cpu_nt.h"
#include "impls/rng_cpu.h"
#include "impls/spqlios_alias.h"
#include "maths_structures.h"
#include "pvda_ffts.h"
#include "rng.h"

#ifdef OPENSSL_BACKEND
#include "rng_openssl.h"
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

void pvda_delete_backend(PvdaBackend* backend)
{
	spqlios_delete_module_info(backend);
	delete_fft_data(backend->pvda_fft_data);
	free(backend);
}

VecUnivDFT* pvda_new_vec_znx_dft(const PvdaBackend* backend, uint64_t size)
{
	return backend->vt.pvda_new_vec_znx_dft(backend, size);
}

int pvda_vec_znx_dft(const PvdaBackend* backend, PolyBivDFT* res, uint64_t res_size, const PolyBiv* a)
{
	return backend->vt.pvda_vec_znx_dft(backend, res, res_size, a);
}

void pvda_delete_vec_znx_dft(const PvdaBackend* backend, double* res)
{
	backend->vt.pvda_delete_vec_znx_dft(backend, res);
}

int64_t* pvda_new_vec_znx_big(const PvdaBackend* backend, int64_t size)
{
	return backend->vt.pvda_new_vec_znx_big(backend, size);
}

int pvda_vec_znx_idft(const PvdaBackend* backend, PolyBiv* res, const PolyBivDFT* a_dft, uint64_t a_size)
{
	return backend->vt.pvda_vec_znx_idft(backend, res, a_dft, a_size);
}

void pvda_delete_vec_znx_big(const PvdaBackend* backend, int64_t* res)
{
	backend->vt.pvda_delete_vec_znx_big(backend, res);
}

PolyUnivDFT* pvda_new_svp_ppol(const PvdaBackend* backend) { return backend->vt.pvda_new_svp_ppol(backend); }

int pvda_svp_prepare(const PvdaBackend* backend, PolyUnivDFT* prepared_pol, const int64_t* pol)
{
	return backend->vt.pvda_svp_prepare(backend, prepared_pol, pol);
}

int pvda_svp_apply_dft(const PvdaBackend* backend, const PolyBivDFT* res, uint64_t res_size,
                       const PolyUnivDFT* prepared_pol, const PolyBiv* a)
{
	return backend->vt.pvda_svp_apply_dft(backend, res, res_size, prepared_pol, a);
}

int pvda_svp_apply_dft_to_dft(const PvdaBackend* backend, const PolyBivDFT* res, uint64_t res_size,
                              const PolyUnivDFT* ppol, const PolyBivDFT* a, uint64_t a_size)
{
	return backend->vt.pvda_svp_apply_dft_to_dft(backend, res, res_size, ppol, a, a_size);
}

void pvda_delete_svp_ppol(const PvdaBackend* backend, PolyUnivDFT* res)
{
	backend->vt.pvda_delete_svp_ppol(backend, res);
}

MatBivDFT* pvda_new_vmp_pmat(const PvdaBackend* backend, uint64_t nrows, uint64_t ncols)
{
	return backend->vt.pvda_new_vmp_pmat(backend, nrows, ncols);
}

int pvda_vmp_prepare_contiguous(const PvdaBackend* backend, MatBivDFT* pmat, const int64_t* mat, uint64_t nrows,
                                uint64_t ncols)
{
	return backend->vt.pvda_vmp_prepare_contiguous(backend, pmat, mat, nrows, ncols);
}

int pvda_vmp_apply_dft(const PvdaBackend* backend, PolyBivDFT* res, uint64_t res_size, const PolyBiv* a,
                       const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols)
{
	return backend->vt.pvda_vmp_apply_dft(backend, res, res_size, a, pmat, nrows, ncols);
}

int pvda_vmp_apply_dft_to_dft(const PvdaBackend* backend, VecBivDFT* res, const uint64_t res_size,
                              const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                              const uint64_t ncols)
{
	return backend->vt.pvda_vmp_apply_dft_to_dft(backend, res, res_size, a_dft, a_size, pmat, nrows, ncols);
}

int pvda_vmp_apply_prepared_to_dft(const PvdaBackend* backend, VecBivDFT* res, const uint64_t res_size,
                                   const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                                   const uint64_t ncols)
{
	return backend->vt.pvda_vmp_apply_prepared_to_dft(backend, res, res_size, a_dft, a_size, pmat, nrows, ncols);
}

void pvda_delete_vmp_pmat(const PvdaBackend* backend, MatBivDFT* pmat)
{
	backend->vt.pvda_delete_vmp_pmat(backend, pmat);
}

int pvda_vec_znx_normalize_base2k(const PvdaBackend* backend, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a)
{
	return backend->vt.pvda_vec_znx_normalize_base2k(backend, log2_base2k, res, a);
}

int pvda_znx_small_product(const PvdaBackend* backend, PolyUniv* res, const PolyUniv* a, const PolyUniv* b)
{
	return backend->vt.pvda_znx_small_product(backend, res, a, b);
}

int pvda_vec_znx_negate(const PvdaBackend* backend, PolyBiv* res, const PolyBiv* a)
{
	return backend->vt.pvda_vec_znx_negate(backend, res, a);
}

int pvda_vec_znx_add(const PvdaBackend* backend, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	return backend->vt.pvda_vec_znx_add(backend, res, a, b);
}

int pvda_vec_znx_sub(const PvdaBackend* backend, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	return backend->vt.pvda_vec_znx_sub(backend, res, a, b);
}

int pvda_znx_automorphism(const PvdaBackend* backend, const int64_t p, PolyUniv* res, const PolyUniv* a)
{
	return backend->vt.pvda_znx_automorphism(backend, p, res, a);
}

int pvda_vec_znx_automorphism(const PvdaBackend* backend, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	return backend->vt.pvda_vec_znx_automorphism(backend, p, res, a);
}

int pvda_vec_znx_rotate(const PvdaBackend* backend, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	return backend->vt.pvda_vec_znx_rotate(backend, p, res, a);
}

uint64_t pvda_module_extract_nn(const PvdaBackend* backend) { return backend->vt.pvda_module_extract_nn(backend); }

int pvda_vmp_prepare_vec(const PvdaBackend* backend, PolyBivPrep* pvec, uint64_t nrows, const PolyBiv* a)
{
	return backend->vt.pvda_vmp_prepare_vec(backend, pvec, nrows, a);
}

int pvda_rand_uniform_pow2(const PvdaBackend* backend, int64_t* result, uint64_t nb_bits)
{
	return backend->vt.pvda_rand_uniform_pow2(backend, result, nb_bits);
}

int pvda_rand_uniform_pow2_vec(const PvdaBackend* backend, int64_t* result, uint64_t n, uint64_t nb_bits)
{
	return backend->vt.pvda_rand_uniform_pow2_vec(backend, result, n, nb_bits);
}

int pvda_rand_uniform_binary_vec(const PvdaBackend* backend, uint64_t* result, uint64_t n)
{
	return backend->vt.pvda_rand_uniform_binary_vec(backend, result, n);
}

int pvda_rand_uniform(const PvdaBackend* backend, int64_t* result, int64_t limit_down, int64_t limit_up)
{
	return backend->vt.pvda_rand_uniform(backend, result, limit_down, limit_up);
}

uint64_t pvda_ring_2nth_root(const PvdaBackend* backend, uint64_t t)
{
	return backend->vt.pvda_ring_2nth_root(backend, t);
}
