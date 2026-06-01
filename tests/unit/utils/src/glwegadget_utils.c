#include "glwegadget_utils.h"
#include <math.h>

int64_t noise_bits_half_prod(const GLWEParams* params_glwe, const GLWEGadgetParams* params_glwegadget) {

		int decomp_noise_bits = params_glwe->kappa * glwe_params_l_b(params_glwe);

		int log2n   = ceil(log2(params_glwe->nn));
		int log2nlt = ceil(log2(params_glwe->nn * params_glwegadget->l_tilde));

		decomp_noise_bits -= log2n;
		decomp_noise_bits -= log2nlt;
		decomp_noise_bits -= 2;

  return decomp_noise_bits;
}
