#ifndef RANDOM_POLYNOMIAL_H
#define RANDOM_POLYNOMIAL_H

#include <stdint.h>

#include "backend.h"
#include "glwe_params.h"
#include "pvda_types.h"

/**
 * @file random_polynomial.h
 *
 * @brief Random polynomials of every prividema polynomial type.
 *
 * Every function draws its randomness from the RNG of the backend (see backend_rng.h), whose implementation is
 * chosen independently of the arithmetic one.
 *
 * Available distributions:
 *
 * | Type          | uniform pow2 | uniform range | binary | ternary | binomial | normal |
 * |---------------|:------------:|:-------------:|:------:|:-------:|:--------:|:------:|
 * | PolyUnivZnx   | x            | x             | x      | x       | x        | x (rounded) |
 * | PolyUnivRnX   |              | x             |        |         |          | x      |
 * | PolyUnivTnX   | x (64 bits)  |               |        |         |          | x      |
 * | PolyUnivZnxMP | x (all bits) |               |        |         |          |        |
 * | PolyBiv       | x            |               |        |         | x        | x      |
 *
 * Univariate functions take the number of coefficients n explicitly, so that they can also fill several
 * consecutive polynomials at once (e.g. the k polynomials of a secret key with n = k * N).
 *
 * All functions return 0 on success and -1 on failure (RNG failure or invalid parameters).
 */

// =============================================
// |                                           |
// |                PolyUnivZnx                |
// |                                           |
// =============================================

/** @brief n coefficients uniform in [-2^(nb_bits-1), 2^(nb_bits-1)), 1 <= nb_bits <= 64. */
int univ_znx_random_uniform_pow2(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, uint64_t nb_bits);

/** @brief n coefficients uniform in [low, high] (both included), low <= high. */
int univ_znx_random_uniform(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, int64_t low, int64_t high);

/** @brief n coefficients uniform in {0, 1}. */
int univ_znx_random_binary(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n);

/** @brief n coefficients uniform in {-1, 0, 1}. */
int univ_znx_random_ternary(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n);

/**
 * @brief n coefficients from the centered binomial distribution of parameter eta: values in [-eta, eta],
 * mean 0, variance eta / 2.
 */
int univ_znx_random_binomial(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, uint64_t eta);

/**
 * @brief n coefficients from the rounded normal distribution: round(x) with x ~ N(0, sigma^2).
 *
 * @retval -1 also if a sample does not fit in an int64_t (sigma too large).
 */
int univ_znx_random_normal(const PvdaBackend* module, PolyUnivZnx* res, uint64_t n, double sigma);

// =============================================
// |                                           |
// |                PolyUnivRnX                |
// |                                           |
// =============================================

/** @brief n coefficients uniform in [low, high), low < high, both finite (53 bits of randomness each). */
int univ_rnx_random_uniform(const PvdaBackend* module, PolyUnivRnX* res, uint64_t n, double low, double high);

/** @brief n coefficients from N(mu, sigma^2). */
int univ_rnx_random_normal(const PvdaBackend* module, PolyUnivRnX* res, uint64_t n, double mu, double sigma);

// =============================================
// |                                           |
// |                PolyUnivTnX                |
// |                                           |
// =============================================

/** @brief n coefficients uniform on the torus (all 64 bits random). */
int univ_tnx_random_uniform(const PvdaBackend* module, PolyUnivTnX* res, uint64_t n);

/**
 * @brief n coefficients from the normal distribution on the torus: x mod 1 with x ~ N(0, sigma^2), rounded to the
 * nearest multiple of 2^-64.
 *
 * @param sigma The standard deviation, on the torus scale (e.g. 2^-20).
 */
int univ_tnx_random_normal(const PvdaBackend* module, PolyUnivTnX* res, uint64_t n, double sigma);

// =============================================
// |                                           |
// |               PolyUnivZnxMP               |
// |                                           |
// =============================================

/** @brief All the 64 * nb_words bits of every coefficient uniformly random. */
int univ_znx_mp_random_uniform(const PvdaBackend* module, PolyUnivZnxMP* res);

// =============================================
// |                                           |
// |                  PolyBiv                  |
// |                                           |
// =============================================

/**
 * @brief Every coefficient of the res->nb_limbs limbs of res uniform in [-2^(nb_bits-1), 2^(nb_bits-1)).
 *
 * Works on strided polynomials (only the nn coefficients of each limb are written).
 */
int biv_random_uniform_pow2(const PvdaBackend* module, PolyBiv* res, uint64_t nb_bits);

/**
 * @brief Computes a random uniform bivariate polynomial.
 *
 * Limbs 0 to precision-1 are sampled uniformly on K bits (normalized digits); the other limbs are not modified.
 *
 * @param module      The backend (its RNG is used).
 * @param params_glwe The bivGLWE parameters (gives K).
 * @param res         The output bivariate polynomial.
 * @param precision   The number of limbs to fill, at most res->nb_limbs.
 */
int biv_random_uniform(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, uint64_t precision);

/**
 * @brief Computes a random normal bivariate polynomial: the decomposition of a \RnX polynomial whose coefficients
 * are sampled from N(0, sigma^2), sigma being the noise standard deviation of the parameters (glwe_params_stdev).
 *
 * @retval -1 also if the noise type of the parameters has no defined standard deviation.
 */
int biv_random_normal(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res);

/**
 * @brief Computes a random binomial bivariate polynomial: the last limb (l_a - 1) is sampled from the centered
 * binomial distribution of parameter eta, the other limbs among the first l_a are set to 0.
 */
int biv_random_binomial(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res, uint64_t eta);

/**
 * @brief Generates size uniform \ZnX polynomials (coefficients in [-2^(nb_bits-1), 2^(nb_bits-1))) directly in the
 * DFT domain.
 */
int biv_dft_random_uniform_pow2(const PvdaBackend* module, BE_PolyBivDFT* res_dft, uint64_t size,
                                uint64_t nb_bits);

// =============================================
// |                                           |
// |                  Noise                    |
// |                                           |
// =============================================

/**
 * @brief Samples the encryption noise of the parameters as a bivariate polynomial.
 *
 * The noise is an integer polynomial e on the last limb of the mask (limb l_a - 1), i.e. a torus noise
 * e * 2^(-K * l_a); every other limb of res is set to 0. Add it to a phase with biv_add. The distribution of e
 * depends on the noise type of the parameters:
 * - NOISE_UNIFORM_POWER_OF_TWO: uniform in [-2^(b-1), 2^(b-1)) with b = fast_uniform_nb_bits;
 * - NOISE_UNIFORM: uniform in [-B, B] with B = uniform_range;
 * - NOISE_BINOMIAL: centered binomial with eta = binomial_n / 2;
 * - NOISE_NORMAL: rounded normal of standard deviation normal_sigma * 2^(K * l_a).
 *
 * @param module      The backend (its RNG is used).
 * @param params_glwe The GLWE parameters (which decide the noise).
 * @param res         The output, with nn coefficients and at least l_a limbs (may be strided).
 */
int biv_random_noise(const PvdaBackend* module, const GLWEParams* params_glwe, PolyBiv* res);

#endif  // RANDOM_POLYNOMIAL_H
