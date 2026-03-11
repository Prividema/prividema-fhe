#include "glwe.h"

#include "logger.h"
#include "rng.h"

//! GLWE PART (begin)
int add_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft)
{
	// GLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);

	// Point to DFT(sk_j * a_j)
	PolyBivDFT* as_j_dft = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j_dft, "as_j_dft's malloc failed in add_mult") < 0) return -1;

	// Point to sk_j * a_j
	PolyBiv* as_j = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j, "as_j's malloc failed in add_mult") < 0) {
		free(as_j_dft);
		return -1;
	}

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++) {
		// The j-ème component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyBiv* a_j               = ct + j * N;

		// Computes DFT(sk_j * a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Computes sk_j * a_j
		if (vec_znx_idft_p(module, as_j, l, as_j_dft, l) < 0)
		{
			free(as_j_dft);
			free(as_j);	
			
			return log_perror("vec_znx_idft_p failed in add_mult");
		}
		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < N * l; p++) {
			res[p] += as_j[p];
		}
	}
	free(as_j_dft);
	free(as_j);

	return 0;
}

int glwe_secret_masking(const MODULE* module, GLWECiphertext* ct, GLWESecretKeyDFT* sk_dft, PolyBiv* phase)
{
	uint64_t N     = ct->params->N;
	uint64_t k     = ct->params->k;
	uint64_t kappa = ct->params->kappa;
	uint64_t l     = poly_biv_size(ct->params);

	if (uniform_random_vec(k * N, ct->vec, l, (k + 1) * N, kappa) < 0)
		return log_perror("uniform_random_vec failed in glwe_secret_masking");

	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	PolyBiv* acc = calloc(N * l, sizeof(int64_t));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_masking.") < 0) return -1;

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	if (add_mult(module, ct->params, acc, ct->vec, sk_dft) < 0) 
	{
		free(acc);
		return -1;
	}

	// Add the phase to acc
	for (uint64_t t = 0; t < N * l; t++) 
		acc[t] += phase[t];

	// The pointer to limb_0(b)
	PolyBiv* b_0 = ct->vec + k * N;

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	if (vec_znx_normalize_base2k_p(module, kappa, b_0, l, (k + 1) * N, acc, l, N) < 0)
	{
		return log_perror("vec_znx_normalize_base2k_p failed in glwe_secret_demasking");
	}
		

	free(acc);

	return 0;
}

int sub_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* res, VecBiv* ct, GLWESecretKeyDFT* sk_dft)
{
	// GLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);

	// Point to DFT(sk_j * a_j)
	PolyBivDFT* as_j_dft = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j_dft, "as_j_dft's malloc failed in sub_mult.") < 0) return -1;

	// Point to sk_j * a_j
	PolyBiv* as_j = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j, "as_j's malloc failed in sub_mult.") < 0) {
		free(as_j_dft);
		return -1;
	}

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++) {
		// The j-ème component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyBiv* a_j               = ct + j * N;

		// Computes DFT(sk_j * a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Computes sk_j * a_j
		if (vec_znx_idft_p(module, as_j, l, as_j_dft, l) < 0)
		{
			free(as_j_dft);
			free(as_j);

			return log_perror("vec_znx_idft_p failed in sub_mult");
		}

		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < N * l; p++) {
			res[p] -= as_j[p];
		}
	}
	free(as_j_dft);
	free(as_j);

	return 0;
}

int glwe_secret_demasking(const MODULE* module, PolyBiv* res, GLWESecretKeyDFT* sk_dft, GLWECiphertext* ct)
{
	// GLWE parameters
	uint64_t N   = ct->params->N;
	uint64_t k   = ct->params->k;
	uint64_t l   = poly_biv_size(ct->params);

	PolyBiv* acc = calloc(N * l, sizeof(int64_t));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_demasking.") < 0) return -1;

	if (sub_mult(module, ct->params, acc, ct->vec, sk_dft) < 0) {
		free(acc);
		log_perror("sub_mult failed in glwe_secret_demasking.");
		return -1;
	}

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	int64_t* b = ct->vec + k * N;
	add_biv_poly(ct->params, acc, N, b, (k + 1) * N, acc, N);

	// The phase in Zn[X,Y]
	if (vec_znx_normalize_base2k_p(module, ct->params->kappa, res, l, N, acc, l, N) < 0)
	{
		free(acc);
		return log_perror("vec_znx_normalize_base2k_p failed in glwe_secret_demasking");
	}

	free(acc);

	return 0;
}

//! GLWE IN DFT SPACE PART (begin)

int glwe_secret_masking_dft(const MODULE* module, GLWECiphertextDFT* ct_dft, GLWESecretKeyDFT* sk_dft, PolyBivDFT* phase_dft)
{
	// GLWE parameters
	const GLWECtParams* params = ct_dft->params;
	uint64_t k           = params->k;
	uint64_t N           = params->N;
	uint64_t kappa       = params->kappa;
	uint64_t l           = poly_biv_size(params);

	VecBiv* ct           = calloc(glwe_coef_number(params), sizeof(int64_t));
	if (log_is_null(ct, "ct's calloc failed in glwe_secret_masking_dft") < 0) return -1;

	if (uniform_random_vec(k * N, ct, l, (k + 1) * N, kappa) < 0) {
		free(ct);
		return -1;
	}

	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	PolyBiv* acc = calloc(N * l, sizeof(int64_t));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_masking_dft.") < 0) {
		free(ct);
		return -1;
	}

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	if (add_mult(module, params, acc, ct, sk_dft) < 0) {
		log_perror("add_mult failed in glwe_secret_masking_dft.");
		free(ct);
		free(acc);
		return -1;
	}

	// The pointer to limb_0(b)
	PolyBiv* b_0 = ct + k * N;

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	if (vec_znx_normalize_base2k_p(module, kappa, b_0, l, (k + 1) * N, acc, l, N) < 0)
		return log_perror("vec_znx_normalize_base2k_p failed in glwe_secret_masking_dft");

	// Computes the GLWE ciphertext in the DFT domain
	vec_znx_dft_p(module, ct_dft->vec, (k + 1) * l, ct, (k + 1) * l, N);

	// Add the phase to the result ciphertext's b
	for (uint64_t i = 1; i <= l; i++) {
		for (uint64_t p = 0; p < N; p++) {
			ct_dft->vec[(i - 1) * (k + 1) * N + k * N + p] += phase_dft[(i - 1) * N + p];
		}
	}

	free(acc);
	free(ct);

	return 0;
}

int glwe_secret_demasking_dft(const MODULE* module, PolyBiv* res, GLWESecretKeyDFT* sk_dft, GLWECiphertextDFT* ct_dft)
{
	// GLWE parameters
	const GLWECtParams* params = ct_dft->params;
	uint64_t N           = params->N;
	uint64_t k           = params->k;
	uint64_t l           = poly_biv_size(params);

	// Computes the input ciphertext out of the DFT domain
	VecBiv* ct = calloc(glwe_coef_number(params), sizeof(int64_t));
	if (log_is_null(ct, "ct's calloc failed in glwe_secret_demasking_dft") < 0) return -1;
	if (vec_znx_idft_p(module, ct, glwe_size(params), ct_dft->vec, glwe_size(params)) < 0)
	{
		free(ct);
		return log_perror("vec_znx_idft_p failed in glwe_secret_demasking_dft");
	}

	PolyBiv* acc = calloc(N * l, sizeof(int64_t));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_demasking_dft.") < 0) {
		free(ct);
		return -1;
	}

	if (sub_mult(module, params, acc, ct, sk_dft) < 0) {
		free(ct);
		free(acc);
		return -1;
	}

	// Computes acc = phi_sk(a,b) = b - Sum_j{0,k-1}[sk_j * a_j]
	PolyBiv* b = ct + k * N;
	add_biv_poly(params, acc, N, b, (k + 1) * N, acc, N);

	// The phase in Zn[X,Y]
	if (vec_znx_normalize_base2k_p(module, params->kappa, res, l, N, acc, l, N) < 0)
	{
		free(acc);
		free(ct);

		return log_perror("vec_znx_normalize_base2k_p failed in glwe_secret_demasking_dft");
	}

	free(acc);
	free(ct);

	return 0;
}
