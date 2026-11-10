#include "backend.h"

#include <stdlib.h>

#include "backend_private.h"
#include "impls/cpu_nt.h"
#include "impls/spqlios_alias.h"
#include "pvda_ffts.h"
#include "logger.h"

PvdaBackend* pvda_new_spqlios_backend(int nn)
{
	PvdaBackend* res = calloc(1, sizeof(PvdaBackend));
	if (!res)
	{
		log_perror("Backend allocation failed");
		return NULL;
	}
	// The RNG is set first: the arithmetic initialization may need randomness
	int rng_status = pvda_backend_set_rng(res, pvda_default_rng_type(), NULL);

	res->spqlios_module = spqlios_new_module_info(nn);
	pvda_fill_spqlios(&res->vt);
	res->pvda_fft_data = new_fft_data(nn);
	pvda_fill_cpu_nt(&res->vt);

	if (rng_status < 0 || !res->spqlios_module || !res->pvda_fft_data)
	{
		log_message(LOG_ERROR, "Backend creation failed");
		pvda_delete_backend(res);
		return NULL;
	}
	return res;
}

void pvda_delete_backend(PvdaBackend* module)
{
	if (!module) return;
	spqlios_delete_module_info(module);
	delete_fft_data(module->pvda_fft_data);
	pvda_delete_backend_rng(module);
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

int pvda_univ_znx_to_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* res, const PolyUnivZnx* a)
{
	return module->vt.pvda_univ_znx_to_dft(module, res, a);
}

int pvda_univ_znx_idft(const PvdaBackend* module, PolyUnivZnx* res, const BE_PolyUnivZnxDFT* pol)
{
	return module->vt.pvda_univ_znx_idft(module, res, pol);
}

// Same pattern as pvda_biv_dft_at: one const vtable slot serves both versions.
BE_PolyBivDFT* pvda_univ_dft_as_biv_dft(const PvdaBackend* module, BE_PolyUnivZnxDFT* pol)
{
	return (BE_PolyBivDFT*)module->vt.pvda_univ_dft_as_biv_dft(module, pol);
}

const BE_PolyBivDFT* pvda_univ_dft_as_biv_dft_const(const PvdaBackend* module, const BE_PolyUnivZnxDFT* pol)
{
	return module->vt.pvda_univ_dft_as_biv_dft(module, pol);
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

int pvda_univ_znx_to_svp_prep(const PvdaBackend* module, BE_PolyUnivSvpPrep* res, const PolyUnivZnx* a)
{
	return module->vt.pvda_univ_znx_to_svp_prep(module, res, a);
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

int pvda_biv_to_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const PolyBiv* a)
{
	return module->vt.pvda_biv_to_dft(module, res, res_size, a);
}

int pvda_biv_idft(const PvdaBackend* module, PolyBiv* res, const BE_PolyBivDFT* pol, uint64_t pol_size)
{
	return module->vt.pvda_biv_idft(module, res, pol, pol_size);
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
                       const BE_PolyUnivSvpPrep* prepared_pol, const PolyBiv* a)
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

int pvda_biv_to_vmp_prep(const PvdaBackend* module, BE_PolyBivVmpPrep* res, uint64_t res_size, const PolyBiv* a)
{
	return module->vt.pvda_biv_to_vmp_prep(module, res, res_size, a);
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

int pvda_biv_to_vmp_pmat(const PvdaBackend* module, BE_PolyMatVmpPrep* pmat, const PolyBiv* mat, uint64_t nrows,
                         uint64_t ncols)
{
	return module->vt.pvda_biv_to_vmp_pmat(module, pmat, mat, nrows, ncols);
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

int pvda_vmp_apply_dft(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_size, const PolyBiv* a,
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

int pvda_biv_normalize_base2k(const PvdaBackend* module, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_biv_normalize_base2k(module, log2_base2k, res, a);
}

int pvda_univ_znx_small_product(const PvdaBackend* module, PolyUnivZnx* res, const PolyUnivZnx* a,
                                const PolyUnivZnx* b)
{
	return module->vt.pvda_univ_znx_small_product(module, res, a, b);
}

int pvda_biv_negate(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_biv_negate(module, res, a);
}

int pvda_biv_add(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	return module->vt.pvda_biv_add(module, res, a, b);
}

int pvda_biv_sub(const PvdaBackend* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	return module->vt.pvda_biv_sub(module, res, a, b);
}

int pvda_univ_znx_automorphism(const PvdaBackend* module, const int64_t p, PolyUnivZnx* res, const PolyUnivZnx* a)
{
	return module->vt.pvda_univ_znx_automorphism(module, p, res, a);
}

int pvda_biv_automorphism(const PvdaBackend* module, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_biv_automorphism(module, p, res, a);
}

int pvda_biv_rotate(const PvdaBackend* module, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_biv_rotate(module, p, res, a);
}

uint64_t pvda_module_extract_nn(const PvdaBackend* module)
{
	return module->vt.pvda_module_extract_nn(module);
}

uint64_t pvda_ring_2nth_root(const PvdaBackend* module, uint64_t t)
{
	return module->vt.pvda_ring_2nth_root(module, t);
}

int pvda_biv_dft_add(const PvdaBackend* module, BE_PolyBivDFT* res, uint64_t res_n_limbs, const BE_PolyBivDFT* a,
                     uint64_t a_n_limbs, const BE_PolyBivDFT* b, uint64_t b_n_limbs)
{
	return module->vt.pvda_biv_dft_add(module, res, res_n_limbs, a, a_n_limbs, b, b_n_limbs);
}