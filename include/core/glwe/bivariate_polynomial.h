#ifndef BIVARIATE_POLYNOMIAL_H
#define BIVARIATE_POLYNOMIAL_H 

#include <stdint.h>
#include "glwe_ct_params.h"
#include "structure_alias.h"

//! BIV POLY PART (begin) 
int64_t poly_biv_coef_number(GLWECtParams* params);
void add_biv_poly(GLWECtParams* params, 
                  PolyBiv* res, int64_t res_sl,
                  PolyBiv* a, int64_t a_sl,
                  PolyBiv* b, int64_t b_sl 
);

//! BIV POLY IN DFT PART (begin)
int64_t poly_biv_coef_number_dft(GLWECtParams* params);
void add_biv_poly_dft(GLWECtParams* params, 
                  PolyBivDFT* res, int64_t res_sl,
                  PolyBivDFT* a, int64_t a_sl,
                  PolyBivDFT* b, int64_t b_sl 
);

//! COMMON PART (begin)
int64_t poly_univ_bytes(GLWECtParams* params);
int64_t poly_biv_size(GLWECtParams* params);
int64_t poly_biv_bytes(GLWECtParams* params);

void biv_to_univ(GLWECtParams* params, double* res_univ, PolyBiv* poly);

#endif // BIVARIATE_POLYNOMIAL_H