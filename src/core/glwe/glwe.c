#include "glwe.h"

#include "glwe_ct_params.h"
#include "logger.h"
#include "rng.h"
#include "utils.h"

//! bivGLWE PART (begin)
int add_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* res, VecBiv* glwe, GLWESecretKeyDFT* sk_dft)
{
	int status = -1;
	// GLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);

	PolyBivDFT* as_j_dft = NULL;  //DFT(sk_j * a_j)
	PolyBiv* as_j        = NULL;  // sk_j * a_j

	as_j_dft = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in add_mult");
	as_j = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j, "as_j's malloc failed in add_mult");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyBiv* a_j               = glwe + j * N;

		// Computes DFT(sk_j * a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Invert DFT to get sk_j * a_j
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l), "vec_znx_idft_p failed in add_mult");

		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < N * l; p++) res[p] += as_j[p];
	}
	status = 0;
cleanup:

	free(as_j_dft);
	free(as_j);

	return status;
}

int sub_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft)
{
	int status = -1;
	// GLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);

	PolyBivDFT* as_j_dft = NULL;  //DFT(sk_j * a_j)
	PolyBiv* as_j        = NULL;  // sk_j * a_j

	as_j_dft = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in sub_mult");
	as_j = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j, "as_j's malloc failed in sub_mult");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyBiv* a_j               = ct + j * N;

		// Computes DFT(sk_j * a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Invert DFT to get sk_j * a_j
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l), "vec_znx_idft_p failed in add_mult");

		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < N * l; p++)
		{
			res[p] -= as_j[p];
		}
	}
	status = 0;
cleanup:

	free(as_j_dft);
	free(as_j);

	return status;
}

int glwe_secret_masking(const MODULE* module, GLWECiphertext* glwe, const GLWESecretKeyDFT* sk_dft,
                        const PolyBiv* phase)
{
	int status = -1;

	const GLWECtParams* params = (const GLWECtParams*)glwe->params;
	VecBiv* result             = glwe->vec;
	// GLWE parameters
	uint64_t N       = params->N;
	uint64_t k       = params->k;
	uint64_t kappa   = params->kappa;
	uint64_t n_limbs = params->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// The acc variable is used for the sum
	// As a reminder, if we deonte acc_j its value at cycle j,
	// at cycle j + 1 we compute
	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	// In other words, acc is <A, SK>
	PolyBiv* acc         = NULL;
	PolyBivDFT* as_j_dft = NULL;  // DFT(sk_j) * DFT(a_j)
	PolyBiv* as_j        = NULL;  // sk_j * a_j

	acc = calloc(N * l, sizeof(double));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking");
	as_j_dft = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_masking");
	as_j = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_masking");

	// Draws uniformly in Zn[X,Y] the ajs'
	CHECK_CALL(uniform_random_vec(k * N, result, l, (k + 1) * N, kappa),
	           "A generation failed in glwe_secret_masking_dft");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of the DFT encoding of the secret key
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];

		// Computes DFT(sk_j) * DFT(a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, result + j * N, l, (k + 1) * N);

		// Undo DFT to retreive sk_j * a_j
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l),
		           "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib");

		// Add it all to the accumulator
		for (uint64_t p = 0; p < N * l; p++) acc[p] += as_j[p];
	}

	// Add the message with noise to acc, so acc now becomes B before normalization
	for (size_t i = 0; i < N * l; ++i) acc[i] += phase[i];

	// The pointer to the last row of the ciphertext vector (B)
	PolyBiv* b_0 = result + k * N;

	// Normalize acc (B) and put it in the result variable
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, b_0, l, (k + 1) * N, acc, l, N),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_masking_ggsw_lib");

	status = 0;

cleanup:
	free(as_j);
	free(as_j_dft);
	free(acc);

	return status;
}

int glwe_secret_demasking(const MODULE* module, PolyBiv* res, const GLWESecretKeyDFT* sk_dft,
                          const GLWECiphertext* glwe)
{
	const GLWECtParams* params = glwe->params;
	const VecBiv* glwe_vec     = glwe->vec;

	int status = -1;

	// GLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);

	// Variables
	PolyBiv* acc         = NULL;  // -Sum_j{0,k-1}[sk_j * a_j]
	PolyBivDFT* as_j_dft = NULL;  // DFT(sk_j * a_j)
	PolyBiv* as_j        = NULL;  // sk_j * a_j

	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_demasking_ggsw_lib");

	as_j_dft = calloc(2 * poly_biv_coef_number_dft(params), sizeof(double));
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_demasking_ggsw_lib");

	as_j = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_demasking_ggsw_lib");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of the secret key in DFT form and the bivGLWE/GLW ciphertext respectively
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		const PolyUniv* a_j        = glwe_vec + j * N;

		// Computes DFT(sk_j * a_j)
		pvda_svp_apply_dft(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Computes sk_j * a_j by inverting the DFT
		CHECK_CALL(pvda_vec_znx_idft(module, as_j, l, as_j_dft, l),
		           "vec_znx_idft_p failed in glwe_secret_demasking_ggsw_lib");

		// And subs it to acc
		for (uint64_t p = 0; p < N * l; p++) acc[p] -= as_j[p];
	}

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	const PolyBiv* b = glwe_vec + k * N;

	// acc += b <=> acc = b - sum(sk_j*a_j)
	add_biv_poly(params, acc, N, acc, N, b, (k + 1) * N);

	//normalize acc into result
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params->kappa, res, l, N, acc, l, N),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_demasking_ggsw_lib");

	status = 0;

cleanup:
	free(as_j);
	free(as_j_dft);
	free(acc);

	return status;
}

//! bivGLWE IN DFT SPACE PART (begin)

int glwe_secret_masking_dft(const MODULE* module, GLWECiphertextDFT* glwe_dft, const GLWESecretKeyDFT* sk_dft,
                            const PolyBivDFT* phase_dft)
{
	int status = -1;
	// GLWE parameters
	const GLWECtParams* params = glwe_dft->params;
	uint64_t k                 = params->k;
	uint64_t N                 = params->N;
	uint64_t kappa             = params->kappa;
	uint64_t l                 = poly_biv_size(params);

	VecBiv* glwe_vec = NULL;
	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	PolyBiv* acc = NULL;
	glwe_vec     = calloc(glwe_coef_number(params), sizeof(int64_t));
	CHECK_ALLOC(glwe_vec, "calloc failed in glwe_secret_masking_dft");
	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "calloc failed in glwe_secret_masking_dft");

	CHECK_CALL(uniform_random_vec(k * N, glwe_vec, l, (k + 1) * N, kappa),
	           "Random vec generation failed in glwe_secret_masking_dft");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	CHECK_CALL(add_mult(module, params, acc, glwe_vec, sk_dft), "add_mult failed in glwe_secret_masking_dft.");

	// The pointer to limb_0(b)
	PolyBiv* b_0 = glwe_vec + k * N;

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, kappa, b_0, l, (k + 1) * N, acc, l, N),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_masking_dft");

	// Computes the bivGLWE ciphertext in the DFT domain
	pvda_vec_znx_dft(module, glwe_dft->vec, (k + 1) * l, glwe_vec, (k + 1) * l, N);

	// Add the phase to the result ciphertext's b
	for (uint64_t i = 0; i < l; i++)
	{
		for (uint64_t p = 0; p < N; p++)
		{
			glwe_dft->vec[i * (k + 1) * N + k * N + p] += phase_dft[i * N + p];
		}
	}

	status = 0;
cleanup:
	free(glwe_vec);
	free(acc);

	return status;
}

int glwe_secret_demasking_dft(const MODULE* module, PolyBiv* res, const GLWESecretKeyDFT* sk_dft,
                              const GLWECiphertextDFT* glwe_dft)
{
	int status = -1;
	// GLWE parameters
	const GLWECtParams* params = glwe_dft->params;
	uint64_t N                 = params->N;
	uint64_t k                 = params->k;
	uint64_t l                 = poly_biv_size(params);

	// Computes the input ciphertext out of the DFT domain
	VecBiv* glwe = NULL;
	PolyBiv* acc = NULL;
	glwe         = calloc(glwe_coef_number(params), sizeof(int64_t));
	CHECK_ALLOC(glwe, "calloc failed in glwe_secret_masking_dft");
	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "calloc failed in glwe_secret_masking_dft");

	CHECK_CALL(pvda_vec_znx_idft(module, glwe, glwe_size(params), glwe_dft->vec, glwe_size(params)),
	           "idft failed in glwe_secret_demasking_dft");

	CHECK_CALL(sub_mult(module, params, acc, glwe, sk_dft), "sub_mult failed in glwe_secret_masking_dft");

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	PolyBiv* b = glwe + k * N;
	add_biv_poly(params, acc, N, b, (k + 1) * N, acc, N);

	// The phase in Zn[X,Y]
	CHECK_CALL(pvda_vec_znx_normalize_base2k(module, params->kappa, res, l, N, acc, l, N),
	           "normalization failed in glwe_secret_demasking_dft");

	status = 0;

cleanup:
	free(acc);
	free(glwe);

	return status;
}
