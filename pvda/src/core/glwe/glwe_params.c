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
		case NOISE_UNIFORM: {
			// Uniform on [-B, B] has variance B(B+1)/3 ~ (B^2)/3
			params_glwe->uniform_range = (int64_t)llround(sqrt(3.0) * rescaled_noise);
			break;
		}
		case NOISE_BINOMIAL: {
			// Centered binomial B(n, 1/2) - n/2 has variance n/4, with n even (n = 2 * eta)
			double n = 4.0 * rescaled_noise * rescaled_noise;
			if (n > 0x1p30) RAISE_ERROR("Noise too big for binomial sampling (n > 2^30), use NOISE_NORMAL instead");
			uint64_t eta            = (uint64_t)ceil(n / 2.0);
			params_glwe->binomial_n = 2 * eta;
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
		case NOISE_UNIFORM: {
			double b = (double)params_glwe->uniform_range;
			return ldexp(sqrt(b * (b + 1.0) / 3.0), -(int)params_glwe->kappa * (int)glwe_params_l_a(params_glwe));
		}
		case NOISE_NORMAL:
			// normal_sigma is stored on the torus scale (see new_glwe_params)
			return params_glwe->normal_sigma;
		case NOISE_BINOMIAL:
			return ldexp(0.5 * sqrt((double)params_glwe->binomial_n),
			             -(int)params_glwe->kappa * (int)glwe_params_l_a(params_glwe));
	}
	return NAN;
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
