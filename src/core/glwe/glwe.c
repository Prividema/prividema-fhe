#include "glwe.h"

#include "logger.h"
#include "utils.h"
#include "rng.h"

//! bivGLWE PART (begin)
int add_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* res, VecBiv* glwe, GLWESecretKeyDFT* sk_dft)
{
	int status = -1;
	
	// bivGLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);
	
	// Variables
	PolyBivDFT* as_j_dft = NULL;
	PolyBiv* as_j = NULL; 
	
	// Point to DFT(sk_j * a_j)
	as_j_dft = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in add_mult");

	// Points to sk_j * a_j
	as_j = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j, "as_j's malloc failed in add_mult");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++) 
	{
		// The j-ème component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyBiv* a_j               = glwe + j * N;

		// Computes DFT(sk_j * a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Computes sk_j * a_j
		CHECK_CALL(vec_znx_idft_p(module, as_j, l, as_j_dft, l), "vec_znx_idft_p failed in add_mult");

		// Computes acc = acc - sk_j * a_j
		for (uint64_t p = 0; p < N * l; p++)
			res[p] += as_j[p];
	}
	
	status = 0;

cleanup:
	free(as_j_dft);
	free(as_j);

	return status;
}

int glwe_secret_masking(const MODULE* module, GLWECiphertext* glwe, GLWESecretKeyDFT* sk_dft, PolyBiv* phase)
{
	int status = -1;

	// Parameters
	uint64_t N     = glwe->params->N;
	uint64_t k     = glwe->params->k;
	uint64_t kappa = glwe->params->kappa;
	uint64_t l     = poly_biv_size(glwe->params);

	// Variables
	PolyBiv* acc = NULL;

	// Draws uniformly in Zn[X] the a_is'
	CHECK_CALL(uniform_random_vec(k * N, glwe->vec, l, (k + 1) * N, kappa) < 0, "uniform_random_vec failed in glwe_secret_masking");

	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking.");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	CHECK_CALL(add_mult(module, glwe->params, acc, glwe->vec, sk_dft), "add_mult failed in glwe_secret_masking");

	// Add the phase to acc
	for (uint64_t t = 0; t < N * l; t++) 
		acc[t] += phase[t];

	// The pointer to limb_0(b)
	PolyBiv* b_0 = glwe->vec + k * N;

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	CHECK_CALL(vec_znx_normalize_base2k_p(module, kappa, b_0, l, (k + 1) * N, acc, l, N), 
			  "vec_znx_normalize_base2k_p failed in glwe_secret_demasking");

	status = 0;

cleanup:
	free(acc);

	return status;
}

int sub_mult(const MODULE* module, const GLWECtParams* params, PolyBiv* res, VecBiv* glwe, GLWESecretKeyDFT* sk_dft)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t N = params->N;
	uint64_t k = params->k;
	uint64_t l = poly_biv_size(params);

	// Variables
	PolyBivDFT* as_j_dft = NULL;
	PolyBiv* as_j = NULL;

	// Point to DFT(sk_j * a_j)
	as_j_dft = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in sub_mult.");

	// Point to sk_j * a_j
	as_j = malloc(poly_biv_bytes(params));
	CHECK_ALLOC(as_j, "as_j's malloc failed in sub_mult.");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++) {
		// The j-ème component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyBiv* a_j               = glwe + j * N;

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
	
	status = 0;

cleanup:
	free(as_j_dft);
	free(as_j);

	return status;
}

int glwe_secret_demasking(const MODULE* module, PolyBiv* res, GLWESecretKeyDFT* sk_dft, GLWECiphertext* glwe)
{
	int status = 0;

	// bivGLWE parameters
	uint64_t N   = glwe->params->N;
	uint64_t k   = glwe->params->k;
	uint64_t l   = poly_biv_size(glwe->params);

	// Variables
	PolyBiv* acc = NULL;

	// Points to phi_sk(a,b) = b - Sum_j{0,k-1}[sk_j * a_j]
	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_demasking.");

	// Computes acc = - Sum_j{0,k-1}[sk_j * a_j]
	CHECK_CALL(sub_mult(module, glwe->params, acc, glwe->vec, sk_dft), "sub_mult failed in glwe_secret_demasking.");

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	int64_t* b = glwe->vec + k * N;
	add_biv_poly(glwe->params, acc, N, b, (k + 1) * N, acc, N);

	// The phase in Zn[X,Y]
	CHECK_CALL(vec_znx_normalize_base2k_p(module, glwe->params->kappa, res, l, N, acc, l, N),
			  "vec_znx_normalize_base2k_p failed in glwe_secret_demasking");
	
	status = 0;

cleanup:
	free(acc);

	return status;
}

//! bivGLWE IN DFT SPACE PART (begin)

int glwe_secret_masking_dft(const MODULE* module, GLWECiphertextDFT* glwe_dft, GLWESecretKeyDFT* sk_dft, PolyBivDFT* phase_dft)
{
	int status = -1;

	// bivGLWE parameters
	const GLWECtParams* params = glwe_dft->params;
	uint64_t k           = params->k;
	uint64_t N           = params->N;
	uint64_t kappa       = params->kappa;
	uint64_t l           = poly_biv_size(params);

	// Parameters
	VecBiv* glwe_vec = NULL;
	PolyBiv* acc = NULL;

	// Points to the glwe_dft's vector out of the DFT domain
	glwe_vec  = calloc(glwe_coef_number(params), sizeof(int64_t));
	CHECK_ALLOC(glwe_vec, "ct's calloc failed in glwe_secret_masking_dft");

	// Draws uniformly in Zn[X] the a_is'
	CHECK_CALL(uniform_random_vec(k * N, glwe_vec, l, (k + 1) * N, kappa), "uniform_random_vec failed in glwe_secret_masking_dft");

	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking_dft.");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	CHECK_CALL(add_mult(module, params, acc, glwe_vec, sk_dft), "add_mult failed in glwe_secret_masking_dft.");

	// The pointer to limb_0(b)
	PolyBiv* b_0 = glwe_vec + k * N;

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	CHECK_CALL(vec_znx_normalize_base2k_p(module, kappa, b_0, l, (k + 1) * N, acc, l, N), 
			  "vec_znx_normalize_base2k_p failed in glwe_secret_masking_dft");

	// Computes the bivGLWE ciphertext in the DFT domain
	vec_znx_dft_p(module, glwe_dft->vec, (k + 1) * l, glwe_vec, (k + 1) * l, N);

	// Add the phase to the result ciphertext's b
	for (uint64_t i = 1; i <= l; i++) {
		for (uint64_t p = 0; p < N; p++) {
			glwe_dft->vec[(i - 1) * (k + 1) * N + k * N + p] += phase_dft[(i - 1) * N + p];
		}
	}

	
	status = 0;

cleanup:
	free(glwe_vec);
	free(acc);

	return status;
}

int glwe_secret_demasking_dft(const MODULE* module, PolyBiv* res, GLWESecretKeyDFT* sk_dft, GLWECiphertextDFT* glwe_dft)
{
	int status = -1;

	// bivGLWE parameters
	const GLWECtParams* params = glwe_dft->params;
	uint64_t N           = params->N;
	uint64_t k           = params->k;
	uint64_t l           = poly_biv_size(params);

	// Variables
	VecBiv* glwe = NULL;
	PolyBiv* acc = NULL;

	// Points to glwe_dft's vec out of the DFT domain
	glwe = calloc(glwe_coef_number(params), sizeof(int64_t));
	CHECK_ALLOC(glwe, "ct's calloc failed in glwe_secret_demasking_dft");
	
	// Computes glwe_dft's vec out of the DFT domain
	CHECK_CALL(vec_znx_idft_p(module, glwe, glwe_size(params), glwe_dft->vec, glwe_size(params)), 
			  "vec_znx_idft_p failed in glwe_secret_demasking_dft");

	// Point to phi_sk(a,b) = b - Sum_j{0,k-1}[sk_j * a_j]
	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_demasking_dft.");
	
	// Computes acc = - Sum_j{0,k-1}[sk_j * a_j]
	CHECK_CALL(sub_mult(module, params, acc, glwe, sk_dft), "sub_mult failed in glwe_secret_demasking_dft");

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	PolyBiv* b = glwe + k * N;
	add_biv_poly(params, acc, N, b, (k + 1) * N, acc, N);

	// The phase in Zn[X,Y]
	CHECK_CALL(vec_znx_normalize_base2k_p(module, params->kappa, res, l, N, acc, l, N), 
			  "vec_znx_normalize_base2k_p failed in glwe_secret_demasking_dft");
	
	status = 0;

cleanup:
	free(acc);
	free(glwe);

	return status;
}
