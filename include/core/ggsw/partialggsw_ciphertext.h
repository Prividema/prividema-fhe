#ifndef PARTIALGGSW_CIPHERTEXT_H
#define PARTIALGGSW_CIPHERTEXT_H

#include "ggsw_ct_params.h"
#include "glwe_ciphertext.h"

typedef struct partialggsw_ciphertext
{
	PartialGGSWCtParams* params;
	int64_t* ct;  // matrix of size n_limbs x l_tilde
} PartialGGSWCiphertext;

// bivGLWE (M/(2^kappa))
// bivGLWE (M/(2^kappa)²)
//...
// bivGLWE (M/(2^kappa)^l_tilde)

#endif  // PARTIALGGSW_CIPHERTEXT_H
