#ifndef UNIVARIATE_POLYNOMIAL_H
#define UNIVARIATE_POLYNOMIAL_H

#include "glwe_params.h"
#include "maths_structures.h"
#include "spqlios_alias.h"

int univ_coefs_to_dft(const MODULE* module, PolyUnivDFT* res_dft, const PolyUniv* a);

int univ_dft_to_coefs(const MODULE* module, PolyUniv* res, const PolyUnivDFT* a_dft);

/**
 * @brief Returns the number of bytes needed to store a univariate polynomial.
 *
 * @param params_glwe The bivGLWE parameters.
 * @return uint64_t
 *
 * @note The number of bytes needed to store an univariate polynomial is the same in and out of the DFT domain.
 */
uint64_t poly_univ_bytes(const GLWEParams* params_glwe);

uint64_t poly_univ_rnx_bytes(const GLWEParams* params_glwe);

PolyUniv* new_univ(const GLWEParams* params_glwe);

PolyUnivDFT* new_univ_dft(const MODULE* module);

PolyUnivRnX* new_univ_rnx(const GLWEParams* params_glwe);

void delete_univ(PolyUniv* pol);

void delete_univ_rnx(PolyUnivRnX* pol);

void delete_univ_dft(PolyUnivDFT* pol);
#endif
