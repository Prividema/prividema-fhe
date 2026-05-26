#include "bivariate_polynomial.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/select.h>
#include <sys/types.h>

#include "glwe_params.h"
#include "logger.h"
#include "rng.h"
#include "univariate_polynomial.h"
#include "utils.h"

//! BIV POLY PART (begin)

uint64_t poly_biv_coef_number(const GLWEParams* params_glwe) { return glwe_params_l_a(params_glwe) * params_glwe->nn; }

PolyBiv* new_biv_poly(const GLWEParams* params_glwe)
{
	PolyBiv* pol = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	CHECK_ALLOC(pol, "pol's malloc failed in new_biv_poly");
	return pol;
cleanup:
	return NULL;
}

int normal_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result)
{
	int status = -1;

	PolyUnivRnX* rd_pol_univ = new_univ_rnx(params_glwe);
	CHECK_ALLOC(rd_pol_univ, "rd_pol_univ's malloc failed.");

	CHECK_CALL(normal_random_vec(rd_pol_univ, params_glwe->nn, 0.0, params_glwe->normal_sigma),
	           "random normal vec generation failed");

	CHECK_CALL(univ_rnx_to_biv(params_glwe, result, rd_pol_univ, 0),
	           "univ to biv conversion failed in normal random biv poly generation");

	status = 0;

cleanup:
	delete_univ_rnx(rd_pol_univ);

	return status;
}

int uniform_random_biv_poly(const GLWEParams* params_glwe, PolyBiv* result, int64_t precision)
{
	int status = -1;

	for (uint64_t i = 1; i <= precision; i++)
		for (uint64_t p = 0; p < params_glwe->nn; p++)
			CHECK_CALL(rand_uniform(result + (i - 1) * params_glwe->nn + p, params_glwe->kappa),
			           "rand_uniform failed in uniform_random_biv_poly.");

	status = 0;

cleanup:

	return status;
}

void add_biv_poly(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a, const PolyBiv* b)
{
	uint64_t nn  = params_glwe->nn;
	uint64_t l_a = glwe_params_l_a(params_glwe);
	pvda_vec_znx_add(module, res, l_a, nn, a, l_a, nn, b, l_a, nn);
}

int add_biv_fast_uni_noise(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a)
{
	int status = -1;

	uint64_t nn = params_glwe->nn;

	int64_t* err = new_univ(params_glwe);

	CHECK_ALLOC(err, "Failed alloc in fast uniform noise generation");

	if (params_glwe->fast_uniform_nb_bits)
		CHECK_CALL(uniform_random_pol_znx(err, nn, params_glwe->fast_uniform_nb_bits),
		           "Uniform noise generation failed");
	else
	{
		log_message(LOG_WARN, "Using error/noise with stdev 0");
		memset(err, 0, nn * sizeof(int64_t));
	}
	PolyUniv* last_limb_res = res + nn * (glwe_params_l_a(params_glwe) - 1);
	PolyUniv* last_limb_a   = a + nn * (glwe_params_l_a(params_glwe) - 1);

	if (res != a) memcpy(res, a, poly_biv_bytes(params_glwe));

	pvda_vec_znx_add(module, last_limb_res, 1, nn, last_limb_a, 1, nn, err, 1, nn);

	status = 0;
cleanup:
	free(err);

	return status;
}

int add_biv_noise(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBiv* a)
{
	switch (params_glwe->noise_type)
	{
		case NOISE_UNIFORM_POWER_OF_TWO:
			return add_biv_fast_uni_noise(module, params_glwe, res, a);
		default:
			RAISE_ERROR("Noise type not implemented");
	}
	return 0;
cleanup:
	return -1;
}

//! BIV POLY IN DFT PART (begin)

uint64_t poly_biv_coef_number_dft(const GLWEParams* params_glwe)
{
	return (glwe_params_l_a(params_glwe) * params_glwe->nn) / 2;
}

PolyBivDFT* new_biv_poly_dft(const GLWEParams* params_glwe)
{
	PolyBivDFT* pol_dft = calloc(poly_biv_coef_number_dft(params_glwe), 2 * sizeof(double));
	CHECK_ALLOC(pol_dft, "pol_dft's malloc failed in new_biv_poly");
	return pol_dft;
cleanup:
	return NULL;
}

//! COMMON PART (begin)

uint64_t poly_biv_bytes(const GLWEParams* params_glwe) { return poly_biv_coef_number(params_glwe) * sizeof(int64_t); }

void biv_to_univ_rnx(const GLWEParams* params_glwe, PolyUnivRnX* res_univ, const PolyBiv* pol_biv)
{
	uint64_t nn      = params_glwe->nn;
	uint64_t kappa   = params_glwe->kappa;
	uint64_t l       = glwe_params_l_a(params_glwe);
	uint64_t start_l = l;

	// res_univ(X^p) = Sum_i{1,l}[poly(X^p, Y^i) * 2^(-kappa*i)]
	double pkappa = exp2(-(double)kappa);
	memset(res_univ, 0, poly_univ_rnx_bytes(params_glwe));
	for (uint64_t i = start_l; i >= 1; --i)
		for (uint64_t p = 0; p < nn; p++)
		{
			res_univ[p] += (double)pol_biv[(i - 1) * nn + p];
			res_univ[p] *= pkappa;
		}
}

int univ_rnx_to_biv_low_precision(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivRnX* pol_univ,
                                  int64_t k_offset)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = glwe_params_l_a(params_glwe);
	uint64_t l_max = INT_ROUND_UP_DIV(53ul, kappa);

	// Fills each pol_biv(X^p, Y^i) with the pol_univ's decomposition coefficients of  in [-2^(kappa* - 1) ; 2^(kappa -
	// 1) - 1]
	int64_t mask = (1LL << kappa) - 1;
	double acc   = 0;
	for (uint64_t i = 1; i <= l && i <= l_max; i++)
	{
		acc += ldexp(1.0, (int)kappa - 1 - (int)kappa * (int)i);
	}

	for (uint64_t p = 0; p < nn; p++)
	{
		double tmp = pol_univ[p] + acc;

		if (tmp < 0) tmp -= floor(tmp);
		for (uint64_t i = 1; i <= k_offset; ++i) res[(i - 1) * nn + p] = 0;
		for (uint64_t i = 1; i + k_offset <= l && i <= l_max; i++)
		{
			res[(i + k_offset - 1) * nn + p] = (((int64_t)ldexp(tmp, i * kappa)) & mask) - (1LL << (kappa - 1));
		}
	}

	status = 0;

cleanup:

	return status;
}

int biv_coefs_to_dft(const MODULE* module, const GLWEParams* params_glwe, PolyBivDFT* res_dft, const PolyBiv* a)
{
	uint64_t nn = params_glwe->nn;
	uint64_t l  = glwe_params_l_a(params_glwe);
	pvda_vec_znx_dft(module, res_dft, l, a, l, nn);
	return 0;
}

int biv_dft_to_coefs(const MODULE* module, const GLWEParams* params_glwe, PolyBiv* res, const PolyBivDFT* a_dft)
{
	uint64_t l = glwe_params_l_a(params_glwe);
	return pvda_vec_znx_idft(module, res, l, a_dft, l);
}

int biv_to_univ_tnx(const GLWEParams* params_glwe, PolyUnivTnX* res_tnx, const PolyBiv* pol)
{
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = glwe_params_l_a(params_glwe);
	uint64_t l_max = INT_ROUND_UP_DIV(64ul, kappa);

	memset(res_tnx, 0, nn * sizeof(*res_tnx));

	for (uint64_t i = 0; i < l && i < l_max; ++i)
	{
		for (uint64_t p = 0; p < nn; p++)
		{
			int shft_amt     = 64 - (int)kappa - (int)(i * kappa);
			uint64_t add_amt = shft_amt > 0 ? ((uint64_t)pol[i * nn + p]) << shft_amt : (pol[i * nn + p] >> -shft_amt);
			res_tnx[p] += add_amt;
		}
	}
	return 0;
}

int univ_tnx_to_biv_low_precision(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivTnX* pol_tnx)
{
	uint64_t nn    = params_glwe->nn;
	uint64_t kappa = params_glwe->kappa;
	uint64_t l     = glwe_params_l_a(params_glwe);
	uint64_t l_max = INT_ROUND_UP_DIV(64ul, kappa);

	uint64_t acc = 0;
	int64_t mask = (1LL << kappa) - 1;
	for (uint64_t i = 1; i <= l && 63 + kappa >= i * kappa; i++)
	{
		acc += 1ULL << (63 + kappa - kappa * i);
	}

	for (uint64_t p = 0; p < nn; p++)
	{
		uint64_t tmp = pol_tnx[p] + acc;

		for (uint64_t i = 1; i <= l && i <= l_max; i++)
		{
			int shft_amt = 64 - (int)(i * kappa);
			if (shft_amt > 0)
				res[(i - 1) * nn + p] = (int64_t)((tmp >> shft_amt) & mask) - (1LL << (kappa - 1));
			else
				res[(i - 1) * nn + p] = (int64_t)((tmp << -shft_amt) & mask) - (1LL << (kappa - 1));
		}
	}

	return 0;
}

/**
 * @brief Helper funciton to fetch bits from a value and put then in the LSB of the return
 *
 * @param num Input value from which to fetch bits
 * @param offset Position of the least significant bit to be fetched. E.g., to fetch the MSB
 *               of num we would have offset 63.
 * @param len Number of bits to fetch
 *
 * @return bits [offset:offset+len] of num in positions [0:len] of the return value, 0 elsewhere
 *
 */
static inline uint64_t select_bits(uint64_t num, uint64_t offset, uint64_t len)
{
	assert(offset <= 63);

	return (num >> offset) & ((1ULL << len) - 1);
}

/**
 * @brief Sign-extend to the full uint64_t from position offset, inverting sign if sgn is 1
 *
 * Intended for "casting" 2-complement values of offset + 1 bits into a full 64-bit word.
 *
 * @param num The input number
 * @param offset The MSB of num from which we want to sign-extend
 * @param sgn Whether to invert the sign of the result (1) or not (0)
 *            Values other than 0 and 1 are undefined
 *
 * @return The result of "casting a offset + 1 bit 2-complement int" to int64_t,
 *         and, if sgn is 1, negating it afterwards
 *
 */
static inline int64_t sgn_ext(uint64_t num, uint64_t offset, int sgn)
{
	uint64_t shift_amt = 63 - offset;
	int64_t s          = ((int64_t)(num << shift_amt)) >> shift_amt;
	return ((s ^ -sgn) + sgn);
}

/**
 * Converts a sign-and-magintude value tnx value (stnx_num) into a bivariate polynomial,
 * putting the LSB of the stnx value in the position of exponent 2^-lsb_pos
 *
 * @param stnx_num A sign-and-magniture tnx number (thus, implicit exponent is 2^-63 and first bit is the sign)
 * @param lsb_pos Position where the LSB of stnx_num should be placed. One should use 63 if no shift is to be performed,
 * but both positive and negative values can be used to multiply/divide by powers of 2 while decomposing
 * @param dst The destination bivariate polynomial
 * @param dst_sl The stride between occurrences of different limbs of the coefficient being decomposed
 * @param params the glwe params (used to get l_a and kappa)
 *
 * @note Left for debugging purposes, as it is used by the rnx conversion function. See biv_decomp_internal_vec
 * for a newer vectorizable (and significantly faster) implementation with better documented internals.
 */
inline static void biv_decomp_internal(uint64_t stnx_num, int lsb_pos, int64_t* dst, int64_t dst_sl,
                                       const GLWEParams* params)
{
	int sgn = (stnx_num & (1UL << 63)) != 0;  //retrieve sign
	stnx_num &= ((1UL << 63) - 1);            //strip sign
	assert(!(stnx_num & (1L << 63)));

	int kappa = (int)params->kappa;
	assert(kappa <= 63);
	int last_l   = INT_ROUND_UP_DIV(lsb_pos, kappa);
	int n_l      = INT_ROUND_UP_DIV(63, kappa);
	int k_offset = last_l * kappa - lsb_pos;
	assert(k_offset >= 0 && k_offset < kappa);

	uint64_t mask = 0;

	for (uint64_t i = 1; i <= n_l && (i * kappa - 1 - k_offset) < 64; ++i)
	{
		mask += 1ULL << (i * kappa - 1 - k_offset);  //hope that the computer optimises this loop
	}

	// Add the mask and xor by it to do the conversion from [0, 2^-K) to [-2^(K-1), -2^(K-1))
	// if we interpret each K bit group as a 2-complement K-long value
	stnx_num += mask;
	stnx_num ^= mask;

	int l_a = glwe_params_l_a(params);

	// zero the limbs that are too insignificant with the lsb_pos provided
	for (int i = l_a - 1; i >= last_l; --i)
	{
		dst[i * dst_sl] = 0;
	}

	int min_i = last_l - l_a < 0 ? 0 : last_l - l_a;
	int maxi1 = last_l - 1;
	int maxi2 = (63 + k_offset) / kappa;
	int max_i = maxi1 < maxi2 ? maxi1 : maxi2;
	for (int i = min_i; i <= max_i; ++i)
	{
		uint64_t s;
		if (i == 0)
		{
			uint64_t stnxt_tmp = stnx_num << k_offset;
			s                  = select_bits(stnxt_tmp, 0, kappa);
		}
		else
		{
			s = select_bits(stnx_num, i * kappa - k_offset, kappa);
		}
		int64_t s2                     = sgn_ext(s, kappa - 1, sgn);
		dst[(last_l - 1 - i) * dst_sl] = s2;
	}
	int min_i_2 = min_i > max_i + 1 ? min_i : max_i + 1;
	for (int i = last_l - 1 - min_i_2; i >= 0; --i)
	{
		dst[i * dst_sl] = 0;
	}
}

// IEEE 754 double precision mask for sign (bit 63) and mantissa (bits[51:0])
#define DOUBLE_SGN_AND_MANTISSA_BMASK ((1UL << 63) + (1UL << 52) - 1)

int univ_rnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivRnX* pol_univ, int64_t bit_offset)
{
	uint64_t nn = params_glwe->nn;
	int kappa   = (int)params_glwe->kappa;
	for (int p = 0; p < nn; ++p)
	{
		double val = pol_univ[p];
		uint64_t s_val;
		// Casts the double bitwise to a uint64_t so we can work
		// directly with the IEEE 754 memory layout
		memcpy(&s_val, &val, sizeof(double));

		// IEEE 754 exponent (select bits and undo 1023 fixed offset)
		int exp = (int)select_bits(s_val, 52, 11);
		exp -= 1023;

		// IEEE 754 significand
		s_val &= DOUBLE_SGN_AND_MANTISSA_BMASK;

		// Add implicit bit (in the case of denormal numbers, this makes the computation be wrong by)
		// the minimum non-denormal distance, which is a somewhat graceful fail behaviour
		s_val |= (1UL << 52);

		biv_decomp_internal(s_val, 52 - exp + bit_offset, res + p, nn, params_glwe);
	}
	return 0;
}

/**
 * Converts a vector of signs and another of mangintudes of tnx values into a bivariate polynomial/vector while
 * putting the LSB of the stnx value in the position of exponent 2^-lsb_pos
 *
 * @param mag_vec A contigious vector of the magnitudes of the input values
 * @param sgn_vec A contigious vector of the signs of the input values
 * @param lsb_pos Position where the LSB of the inputs should be placed. One should use 63 if no shift is to be performed,
 * but both positive and negative values can be used to multiply/divide by powers of 2 while decomposing
 * @param dst The destination bivariate polynomial
 * @param dst_sl The stride between occurrences of different limbs of the same coefficient in the output.
 * @param params the glwe params (used to get l_a and kappa)
 *
 *
 */
inline static void biv_decomp_internal_vec(uint64_t* mag_vec, uint8_t* sgn_vec, int lsb_pos, int64_t* dst,
                                           int64_t dst_sl, const GLWEParams* params)
{
	uint64_t nn = params->nn;

	int kappa = (int)params->kappa;
	assert(kappa <= 63);

	int last_l   = INT_ROUND_UP_DIV(lsb_pos, kappa);  // limb number wehre lsb_pos falls
	int k_offset = last_l * kappa - lsb_pos;          // position (bit) in the last written limb where lsb_pos falls
	// Equivalently, number of bits in the last limb that are to the right of it and thus will be 0
	assert(k_offset >= 0 && k_offset < kappa);

	uint64_t mask = 0;

	// Compute the offset-carry mask that converts [0, 2^-K) into [-2^(K-1), -2^(K-1))
	for (uint64_t i = 1; (i * kappa - 1 - k_offset) < 64; ++i)
	{
		mask += 1ULL << (i * kappa - 1 - k_offset);
	}

	// Add the mask and xor by it to do the conversion from [0, 2^-K) to [-2^(K-1), -2^(K-1))
	// if we interpret each K bit group as a 2-complement K-long value
	//
	// Later, the sgn_ext function will be used to retrieve each K-bit group into its own
	// bivariate coefficient, negating it if its sign is negative
	for (int p = 0; p < nn; ++p)
	{
		mag_vec[p] += mask;
		mag_vec[p] ^= mask;
	}

	int l_a = glwe_params_l_a(params);

	// We should start at either limb last_l (last_l-1 position)
	// or l_a if l_a is smaller
	int min_i = last_l - l_a < 0 ? 0 : last_l - l_a;  // min_i = max (last_l - l_a, 0)

	// We should stop at either:
	int maxi1 = last_l - 1;               // Last limb corresponds to position 0, limb 1, most significant limb
	int maxi2 = (63 + k_offset) / kappa;  // Number of limbs that the magnitude can have non-zero values for,
	                                      // as retrieved from the number of bits used
	int max_i = maxi1 < maxi2 ? maxi1 : maxi2;
	for (int i = min_i; i <= max_i; ++i)
	{
		if (i == 0)
			for (int p = 0; p < nn; ++p)
			{
				// For the least significant limb, take into account the rightmost bits that need to be set to 0 due to
				// lsb_pos/k_offset
				uint64_t stnxt_tmp                 = mag_vec[p] << k_offset;
				uint64_t s                         = select_bits(stnxt_tmp, 0, kappa);
				int64_t s2                         = sgn_ext(s, kappa - 1, sgn_vec[p]);
				dst[(last_l - 1 - i) * dst_sl + p] = s2;
			}
		else
			for (int p = 0; p < nn; ++p)
			{
				// Select the correspoding K bit group
				uint64_t s = select_bits(mag_vec[p], i * kappa - k_offset, kappa);
				// Extend the K-bit 2-complement value (that the mask generated) into 64 bits,
				// and invert the sign if the input sign was negative
				// This step is what makes the output be in [-2^(K-1), 2^(K-1)] instead of [-2^(K-1), 2^(K-1))
				int64_t s2                         = sgn_ext(s, kappa - 1, sgn_vec[p]);
				dst[(last_l - 1 - i) * dst_sl + p] = s2;
			}
	}
}

int univ_tnx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUnivTnX* pol_univ, int64_t bit_offset)
{
	int status  = -1;
	uint64_t nn = params_glwe->nn;
	int kappa   = (int)params_glwe->kappa;

	uint64_t* mag_vec = malloc(sizeof(uint64_t) * nn);
	uint8_t* sgn_vec  = malloc(sizeof(uint8_t) * nn);
	CHECK_ALLOC(mag_vec, "Failed malloc in tnx biv conversion");
	CHECK_ALLOC(sgn_vec, "Failed malloc in tnx biv conversion");

	memset(res, 0, poly_biv_bytes(params_glwe));

	// Decompose the tnx values into sign and magintude
	for (int p = 0; p < nn; ++p)
	{
		uint64_t tnx_val = pol_univ[p];
		// The mask is all 1s for negative numbers, all 0 for positive
		uint64_t mask = ((int64_t)tnx_val) >> 63;
		// Branchless absolute value using the mask
		uint64_t abs_val = (tnx_val ^ mask) - mask;
		// magnitude (divided by 2 to ensure no carry overflow)
		mag_vec[p] = (abs_val >> 1);
		sgn_vec[p] = mask & 1;
	}
	// lsb_pos is 63 + bit_offset due to mag_vec being a tnx to which we have removed its LSB
	// and shifted right by 1, thus its implicit denominator is 2^-63
	biv_decomp_internal_vec(mag_vec, sgn_vec, 63 + bit_offset, res, nn, params_glwe);
	status = 0;
cleanup:
	free(mag_vec);
	free(sgn_vec);
	return status;
}

int univ_znx_to_biv(const GLWEParams* params_glwe, PolyBiv* res, const PolyUniv* pol_univ, int64_t bit_offset)
{
	int status = -1;
	if (bit_offset <= 0)
	{
		log_message(LOG_INFO,
		            "A ZnX polynomial is being converted to bivariate with negative or 0 bit_offset. "
		            "The result will always be 0 in that case");
	}
	uint64_t nn = params_glwe->nn;
	int kappa   = (int)params_glwe->kappa;

	uint64_t* mag_vec = malloc(sizeof(uint64_t) * nn);
	uint8_t* sgn_vec  = malloc(sizeof(uint8_t) * nn);
	CHECK_ALLOC(mag_vec, "Failed malloc in tnx biv conversion");
	CHECK_ALLOC(sgn_vec, "Failed malloc in tnx biv conversion");

	memset(res, 0, poly_biv_bytes(params_glwe));

	// Decompose the znx values into sign and magintude
	for (int p = 0; p < nn; ++p)
	{
		uint64_t tnx_val = pol_univ[p];
		// The mask is all 1s for negative numbers, all 0 for positive
		uint64_t mask = ((int64_t)tnx_val) >> 63;
		// Branchless absolute value using the mask
		uint64_t abs_val = (tnx_val ^ mask) - mask;
		// magnitude, not divided since we are in znx and not tnx.
		// The user should make sure that the proper bounds are used.
		mag_vec[p] = abs_val;
		sgn_vec[p] = mask & 1;
	}
	biv_decomp_internal_vec(mag_vec, sgn_vec, bit_offset, res, nn, params_glwe);
	status = 0;
cleanup:
	free(mag_vec);
	free(sgn_vec);
	return status;
}
