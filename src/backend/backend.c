#include "backend.h"

#include <stdlib.h>

#include "backend_private.h"
#include "spqlios_alias.h"

void pvda_fill_spqlios(struct pvda_virtual_table* vt)
{
	vt->pvda_new_vec_znx_dft           = spqlios_new_vec_znx_dft;
	vt->pvda_vec_znx_dft               = spqlios_vec_znx_dft;
	vt->pvda_delete_vec_znx_dft        = spqlios_delete_vec_znx_dft;
	vt->pvda_new_vec_znx_big           = spqlios_new_vec_znx_big;
	vt->pvda_vec_znx_idft              = spqlios_vec_znx_idft;
	vt->pvda_delete_vec_znx_big        = spqlios_delete_vec_znx_big;
	vt->pvda_new_svp_ppol              = spqlios_new_svp_ppol;
	vt->pvda_svp_prepare               = spqlios_svp_prepare;
	vt->pvda_svp_apply_dft             = spqlios_svp_apply_dft;
	vt->pvda_svp_apply_dft_to_dft      = spqlios_svp_apply_dft_to_dft;
	vt->pvda_delete_svp_ppol           = spqlios_delete_svp_ppol;
	vt->pvda_new_vmp_pmat              = spqlios_new_vmp_pmat;
	vt->pvda_vmp_prepare_contiguous    = spqlios_vmp_prepare_contiguous;
	vt->pvda_vmp_apply_dft             = spqlios_vmp_apply_dft;
	vt->pvda_vmp_apply_dft_to_dft      = spqlios_vmp_apply_dft_to_dft;
	vt->pvda_vmp_apply_prepared_to_dft = spqlios_vmp_apply_prepared_to_dft;
	vt->pvda_delete_vmp_pmat           = spqlios_delete_vmp_pmat;
	vt->pvda_vec_znx_normalize_base2k  = spqlios_vec_znx_normalize_base2k;
	vt->pvda_znx_small_product         = spqlios_znx_small_product;
	vt->pvda_vec_znx_negate            = spqlios_vec_znx_negate;
	vt->pvda_vec_znx_add               = spqlios_vec_znx_add;
	vt->pvda_vec_znx_sub               = spqlios_vec_znx_sub;
	vt->pvda_znx_automorphism          = spqlios_znx_automorphism;
	vt->pvda_vec_znx_automorphism      = spqlios_vec_znx_automorphism;
	vt->pvda_vec_znx_rotate            = spqlios_vec_znx_rotate;
	vt->pvda_module_extract_nn         = spqlios_module_extract_nn;
	vt->pvda_vmp_prepare_vec           = spqlios_vmp_prepare_vec;
}

PvdaModule* pvda_new_module_info(int nn)
{
	PvdaModule* res     = malloc(sizeof(PvdaModule));
	res->spqlios_module = spqlios_new_module_info(nn);
	pvda_fill_spqlios(&res->vt);
	return res;
}

void pvda_delete_module_info(PVDA_MODULE* module)
{
	//TODO
}

VecUnivDFT* pvda_new_vec_znx_dft(const PvdaModule* module, uint64_t size)
{
	return module->vt.pvda_new_vec_znx_dft(module, size);
}

int pvda_vec_znx_dft(const PvdaModule* module, double* res, uint64_t res_size, const PolyBiv* a)
{
	return module->vt.pvda_vec_znx_dft(module, res, res_size, a);
}

void pvda_delete_vec_znx_dft(const PvdaModule* module, double* res) { module->vt.pvda_delete_vec_znx_dft(module, res); }

int64_t* pvda_new_vec_znx_big(const PvdaModule* module, int64_t size)
{
	return module->vt.pvda_new_vec_znx_big(module, size);
}

int pvda_vec_znx_idft(const PvdaModule* module, PolyBiv* res, const double* a_dft, uint64_t a_size)
{
	return module->vt.pvda_vec_znx_idft(module, res, a_dft, a_size);
}

void pvda_delete_vec_znx_big(const PvdaModule* module, int64_t* res)
{
	module->vt.pvda_delete_vec_znx_big(module, res);
}

double* pvda_new_svp_ppol(const PvdaModule* module) { return module->vt.pvda_new_svp_ppol(module); }

int pvda_svp_prepare(const PvdaModule* module, PolyUnivDFT* prepared_pol, const int64_t* pol)
{
	return module->vt.pvda_svp_prepare(module, prepared_pol, pol);
}

int pvda_svp_apply_dft(const PvdaModule* module, const double* res, uint64_t res_size, const PolyUnivDFT* prepared_pol,
                       const PolyBiv* a)
{
	return module->vt.pvda_svp_apply_dft(module, res, res_size, prepared_pol, a);
}

int pvda_svp_apply_dft_to_dft(const PvdaModule* module, const double* res, uint64_t res_size, const PolyUnivDFT* ppol,
                              const PolyBivDFT* a, uint64_t a_size)
{
	return module->vt.pvda_svp_apply_dft_to_dft(module, res, res_size, ppol, a, a_size);
}

void pvda_delete_svp_ppol(const PvdaModule* module, double* res) { module->vt.pvda_delete_svp_ppol(module, res); }

double* pvda_new_vmp_pmat(const PvdaModule* module, uint64_t nrows, uint64_t ncols)
{
	return module->vt.pvda_new_vmp_pmat(module, nrows, ncols);
}

int pvda_vmp_prepare_contiguous(const PvdaModule* module, double* pmat, const int64_t* mat, uint64_t nrows,
                                uint64_t ncols)
{
	return module->vt.pvda_vmp_prepare_contiguous(module, pmat, mat, nrows, ncols);
}

int pvda_vmp_apply_dft(const PvdaModule* module, double* res, uint64_t res_size, const PolyBiv* a,
                       const MatBivDFT* pmat, uint64_t nrows, uint64_t ncols)
{
	return module->vt.pvda_vmp_apply_dft(module, res, res_size, a, pmat, nrows, ncols);
}

int pvda_vmp_apply_dft_to_dft(const PvdaModule* module, VecBivDFT* res, const uint64_t res_size, const VecBivDFT* a_dft,
                              uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows, const uint64_t ncols)
{
	return module->vt.pvda_vmp_apply_dft_to_dft(module, res, res_size, a_dft, a_size, pmat, nrows, ncols);
}

int pvda_vmp_apply_prepared_to_dft(const PvdaModule* module, VecBivDFT* res, const uint64_t res_size,
                                   const VecBivDFT* a_dft, uint64_t a_size, const MatBivDFT* pmat, const uint64_t nrows,
                                   const uint64_t ncols)
{
	return module->vt.pvda_vmp_apply_prepared_to_dft(module, res, res_size, a_dft, a_size, pmat, nrows, ncols);
}

void pvda_delete_vmp_pmat(const PvdaModule* module, double* pmat) { module->vt.pvda_delete_vmp_pmat(module, pmat); }

int pvda_vec_znx_normalize_base2k(const PvdaModule* module, uint64_t log2_base2k, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_vec_znx_normalize_base2k(module, log2_base2k, res, a);
}

int pvda_znx_small_product(const PvdaModule* module, PolyUniv* res, const PolyUniv* a, const PolyUniv* b)
{
	return module->vt.pvda_znx_small_product(module, res, a, b);
}

int pvda_vec_znx_negate(const PvdaModule* module, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_vec_znx_negate(module, res, a);
}

int pvda_vec_znx_add(const PvdaModule* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	return module->vt.pvda_vec_znx_add(module, res, a, b);
}

int pvda_vec_znx_sub(const PvdaModule* module, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	return module->vt.pvda_vec_znx_sub(module, res, a, b);
}

int pvda_znx_automorphism(const PvdaModule* module, const int64_t p, PolyUniv* res, const PolyUniv* a)
{
	return module->vt.pvda_znx_automorphism(module, p, res, a);
}

int pvda_vec_znx_automorphism(const PvdaModule* module, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_vec_znx_automorphism(module, p, res, a);
}

int pvda_vec_znx_rotate(const PvdaModule* module, const int64_t p, PolyBiv* res, const PolyBiv* a)
{
	return module->vt.pvda_vec_znx_rotate(module, p, res, a);
}

uint64_t pvda_module_extract_nn(const PvdaModule* module) { return module->vt.pvda_module_extract_nn(module); }

int pvda_vmp_prepare_vec(const PvdaModule* module, double* pvec, uint64_t nrows, const PolyBiv* a)
{
	return module->vt.pvda_vmp_prepare_vec(module, pvec, nrows, a);
}
