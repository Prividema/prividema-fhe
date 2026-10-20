#include "glwe_params.h"

#include <math.h>
#include <stdlib.h>

#include "utils.h"

GLWEParams* new_glwe_params(uint64_t nn, uint64_t k, uint64_t kappa, uint64_t nb_limbs, double sigma,
                            NoiseType noise_type)
{
	GLWEParams* params_glwe = malloc(sizeof(GLWEParams));
	CHECK_ALLOC(params_glwe, "params' malloc failed in new_glwe_ct_params");

	params_glwe->nn                  = nn;
	params_glwe->k                   = k;
	params_glwe->kappa               = kappa;
	params_glwe->ciphertext_nb_limbs = nb_limbs;
	params_glwe->noise_type          = noise_type;

	double rescaled_noise = ldexp(sigma, (int)kappa * (int)glwe_params_l_a(params_glwe));

	if (rescaled_noise > ldexp(1.0, 63)) RAISE_ERROR("Noise did not fit in the bivariate parameters (too big)");
	if (rescaled_noise < 1 && sigma != 0) RAISE_ERROR("Noise did not fit in bivariate parameters (too small)");

	switch (noise_type)
	{
		case NOISE_UNIFORM_POWER_OF_TWO: {
			// Placeholder sigma
			if (sigma == 0)
			{
				params_glwe->fast_uniform_nb_bits = 0;
				break;
			}
			double range                 = sqrt(3.0) * rescaled_noise;
			uint64_t log_2               = ceil(log2(range));
			params_glwe->fast_uniform_nb_bits = log_2;
			break;
		}
		case NOISE_NORMAL: {
			params_glwe->normal_sigma = sigma;
			break;
		}
		default:
			RAISE_ERROR("Trying to use a non-implemented noise type");
	}

	return params_glwe;
cleanup:
	free(params_glwe);
	return NULL;
}

void delete_glwe_params(GLWEParams* params_glwe) { free(params_glwe); }

uint64_t glwe_params_l_a(const GLWEParams* params_glwe)
{
	return (params_glwe->ciphertext_nb_limbs + 1) / (params_glwe->k + 1);
}

uint64_t glwe_params_l_b(const GLWEParams* params_glwe)
{
	return params_glwe->ciphertext_nb_limbs - params_glwe->k * glwe_params_l_a(params_glwe);
}

uint64_t glwe_params_n_limbs(const GLWEParams* params_glwe) { return params_glwe->ciphertext_nb_limbs; }

double glwe_params_stdev(const GLWEParams* params_glwe)
{
	switch (params_glwe->noise_type)
	{
		case NOISE_UNIFORM_POWER_OF_TWO:
			return NAN;
		case NOISE_UNIFORM:
			return NAN;
		case NOISE_NORMAL:
			return ldexp(params_glwe->normal_sigma, -(int)params_glwe->kappa * (int)glwe_params_l_a(params_glwe));
		case NOISE_BINOMIAL:
			return ldexp(0.5 * sqrt(params_glwe->binomial_n),
			             -(int)params_glwe->kappa * (int)glwe_params_l_a(params_glwe));
	}
}

uint64_t glwe_params_bytes(const GLWEParams* params_glwe)
{
	uint64_t nn = params_glwe->nn;
	return glwe_params_n_limbs(params_glwe) * nn * sizeof(int64_t);
}
uint64_t glwe_coef_number(const GLWEParams* params_glwe) { return glwe_params_n_limbs(params_glwe) * params_glwe->nn; }

uint64_t glwe_coef_number_dft(const GLWEParams* params_glwe)
{
	return glwe_params_n_limbs(params_glwe) * params_glwe->nn / 2;
}

double glwe_bivariate_epsilon(const GLWEParams* params_glwe)
{
	return ldexp(1.0, -glwe_params_l_b(params_glwe) * params_glwe->kappa);
}
