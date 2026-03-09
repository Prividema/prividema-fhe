#include "ggsw.h"

#include <stdint.h>
#include <string.h>

#include "glwe_ciphertext.h"
#include "logger.h"
#include "math.h"
#include "rng.h"
#include "spqlios_alias.h"
#include "utils.h"

//! bivGGSW PART (begin)

int add_bivariate_error(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result, const PolyBiv* pol)
{
	int status = -1;

	// Variables 
	PolyBiv* err = NULL;

	// Draw a random error
	err = new_biv_poly(params_glwe);
    CHECK_ALLOC(err, "new_biv_poly failed in add_error");
	CHECK_CALL(normal_random_biv_poly(params_glwe, err), "normal_random_biv_poly failed in add_error");

	// Add the error in the DFT domain
	add_biv_poly(params_glwe, result, params_glwe->N, pol, params_glwe->N, err, params_glwe->N);

	status = 0;

cleanup:
	free(err);

	return status;
}

int glwe_secret_demasking_ggsw_lib(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result,
                                   const GGSWSecretKeyDFT* sk_dft, const VecBiv* glwe_vec)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t N = params_glwe->N;
	uint64_t k = params_glwe->k;
	uint64_t l = poly_biv_size(params_glwe);

	// Variables
	PolyBiv* acc         = NULL;  // -Sum_j{0,k-1}[sk_j * a_j]
	PolyBivDFT* as_j_dft = NULL;  // DFT(sk_j * a_j)
	PolyBiv* as_j        = NULL;  // sk_j * a_j

	acc = calloc(N * l, sizeof(int64_t));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_demasking_ggsw_lib");

	as_j_dft = calloc(2 * poly_biv_coef_number_dft(params_glwe), sizeof(double));
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_demasking_ggsw_lib");

	as_j = calloc(poly_biv_coef_number(params_glwe), sizeof(int64_t));
	CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_demasking_ggsw_lib");

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-ème component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		const PolyUniv* a_j        = glwe_vec + j * N;

		// Computes DFT(sk_j * a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Computes sk_j * a_j by inverting the DFT
		CHECK_CALL(vec_znx_idft_p(module, as_j, l, as_j_dft, l),
		           "vec_znx_idft_p failed in glwe_secret_demasking_ggsw_lib");

		// And subs it to acc
		for (uint64_t p = 0; p < N * l; p++) acc[p] -= as_j[p];
	}

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	const PolyBiv* b = glwe_vec + k * N;

	// acc += b <=> acc = b - sum(sk_j*a_j)
	add_biv_poly(params_glwe, acc, N, acc, N, b, (k + 1) * N);

	//normalize acc into result
	CHECK_CALL(vec_znx_normalize_base2k_p(module, params_glwe->kappa, result, l, N, acc, l, N), 
		"vec_znx_normalize_base2k_p failed in glwe_secret_demasking_ggsw_lib");

	status = 0;

cleanup:
	free(as_j);
	free(as_j_dft);
	free(acc);

	return status;
}

int glwe_secret_masking_ggsw_lib(const MODULE* module, const GLWECtParams* params_glwe, VecBiv* result,
                                 const GGSWSecretKeyDFT* sk_dft, const PolyBiv* phase)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t kappa   = params_glwe->kappa;
	uint64_t n_limbs = params_glwe->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// Variables
	PolyBiv* acc         = NULL;
	PolyBivDFT* as_j_dft = NULL;  // DFT(sk_j) * DFT(a_j)
	PolyBiv* as_j        = NULL;

	// Draws uniformly in Zn[X,Y] the ajs'
	if (uniform_random_vec(k * N, result, l, (k + 1) * N, kappa) < 0) return -1;

	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	// In other words, acc is <A, SK>
	acc = calloc(N * l, sizeof(double));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking_ggsw_lib");

	as_j_dft = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(as_j_dft, "as_j_dft's calloc failed in glwe_secret_masking_ggsw_lib");

	// Point to sk_j * a_j
	// Allocate the variable for sk_j * a_j
	as_j = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(as_j, "as_j's calloc failed in glwe_secret_masking_ggsw_lib");

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-th component of the DFT encoding of the secret key
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];

		// Computes DFT(sk_j) * DFT(a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, result + j * N, l, (k + 1) * N);

		// Undo DFT to retreive sk_j * a_j
		CHECK_CALL(vec_znx_idft_p(module, as_j, l, as_j_dft, l),
		           "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib");

		// Add it all to the accumulator
		for (uint64_t p = 0; p < N * l; p++) acc[p] += as_j[p];
	}

	// Add the message with noise to acc, so acc now becomes B before normalization
	for (size_t i = 0; i < N * l; ++i) acc[i] += phase[i];

	// The pointer to the last row of the ciphertext vector (B)
	PolyBiv* b_0 = result + k * N;

	// Normalize acc (B) and put it in the result variable
	CHECK_CALL(vec_znx_normalize_base2k_p(module, kappa, b_0, l, (k + 1) * N, acc, l, N),
	           "vec_znx_normalize_base2k_p failed in glwe_secret_masking_ggsw_lib");

	status = 0;

cleanup:
	free(as_j);
	free(as_j_dft);
	free(acc);

	return status;
}


int ggsw_secret_encrypt(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertext* result,
                        const GGSWSecretKeyDFT* sk_dft, const PolyUniv* m_univ)
{
	if (params_ggsw->k_tilde > params_ggsw->params_glwe->k)
		return log_perror("k_tilde should not be greater than k in ggsw_secret_encrypt");

	int status = -1;

	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	// bivGLWE parameters
	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t k_tilde = params_ggsw->k_tilde;

	// Variables
	PolyUnivDFT* m_univ_dft   = NULL;  // DFT(msg)
	PolyUnivDFT* m_skj_univ_dft = NULL;  // DFT(msg * sk_j)
	PolyUniv* m_skj_univ        = NULL;  // -msg * sk_j

	// compute_phase_ij requires some extra temp space
	// TODO: find out why it is being set to poly_univ_bytes even though it is a different data type
	double* tmp_sp1 = NULL;

	// Temp space for -m * sk * 2^{-kappa_tilde}
	PolyBiv* glwe_biv_msg = NULL;

	m_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_univ_dft, "malloc failed in ggsw_secret_encrypt");
	m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ_dft, "malloc failed in ggsw_secret_encrypt");
	m_skj_univ = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ, "malloc failed in ggsw_secret_encrypt");
	tmp_sp1 = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(tmp_sp1, "malloc failed in ggsw_secret_encrypt");
	glwe_biv_msg = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg, "malloc failed in ggsw_secret_encrypt");

	// Computes DFT(msg)
	vec_znx_dft_p(module, m_univ_dft, 1, m_univ, 1, N);

	for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
	{
		for (uint64_t j = 0; j < k_tilde + 1; j++)
		{
			{
				// For j < k, m_skj_univ is -m * skj
				// For j = k  we skip this and we get it from m_univ directly
				if (j < params_glwe->k)
				{
					// Computes DFT(msg * sk_j)
					mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, m_univ_dft, 1);

					// Computes -DFT(msg * sk_j)
					for (uint64_t p = 0; p < N; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];

					// Invert the DFT to get -msg * sk_j
					CHECK_CALL(vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1),
					           "vec_znx_idft_p failed in compute_phase_ij");
				}

				// Computes m_skj_univ / 2^{kappa_tilde*i}
				for (uint64_t p = 0; p < N; p++)
					tmp_sp1[p] =
					    ldexp((k == j) ? (double)m_univ[p] : (double)m_skj_univ[p], -params_ggsw->kappa_tilde * i);

				// Compute the base-2^kappa decomposition of tmp_sp1
				CHECK_CALL(univ_to_biv(params_glwe, glwe_biv_msg, tmp_sp1), "univ_to_biv failed in compute_phase_ij");

				// Computes Dec_Kappa(val / 2^{kappa_tilde*i}) + err
				CHECK_CALL(add_bivariate_error(module, params_glwe, glwe_biv_msg, glwe_biv_msg),
				           "add_error failed in compute_phase_ij");
			}
			// Get the pointer for the result position
			VecBiv* glwe_vec = ggsw_retreive_bivglwe(params_ggsw, result->mat, j, i);

			//Compute: bivGLWE(glwe_biv_msg) into glwe_vec
			CHECK_CALL(glwe_secret_masking_ggsw_lib(module, params_glwe, glwe_vec, sk_dft, glwe_biv_msg),
			           "glwe_secret_masking_ggsw_lib failed in ggsw_secret_encrypt");
		}
	}

	status = 0;

cleanup:
	free(glwe_biv_msg);
	free(tmp_sp1);
	free(m_skj_univ_dft);
	free(m_skj_univ);
	free(m_univ_dft);

	return status;
}

int ggsw_external_product(const MODULE* module,
                          GLWECiphertext* result,      // result
                          const GLWECiphertext* glwe,  // bivGLWE ciphertext
                          const GGSWCiphertext* ggsw   // bivGGSW ciphertext
)
{
	int status = -1;

	// Degree of chosen cyclotomic polynomial
	uint64_t N = result->params->N;

	// The bivGGSW ciphertext ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// The bivGLWE ciphertext glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// As the result of the vector-matrix product glwe * ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ggsw->params->n_limbs_tilde;
	uint64_t ncols = ggsw->params->params_glwe->n_limbs;

	// Variables 
	MatBivDFT* ggsw_pmat = NULL;
	VecBivDFT* result_dft = NULL;

	// Point to the bivGGSW ciphertext prepared in the DFT domain
	ggsw_pmat = malloc(ggsw_bytes(ggsw->params));
	CHECK_ALLOC(ggsw_pmat, "mat_dft's malloc failed in ggsw_external_product");

	// Prepares bivGGSW ciphertext prepared in the DFT domain 
	CHECK_CALL(vmp_prepare_contiguous_p(module, ggsw_pmat, ggsw->mat, nrows, ncols), 
		"vmp_prepare_contiguous_p failed in ggsw_external_product");

	// The pointer to ExternalProduct(glwe, ggsw)
	result_dft = malloc(glwe_bytes(ggsw->params->params_glwe));
	CHECK_ALLOC(result_dft, "result's malloc failed in ggsw_external_product");

	// Computes ExternalProduct(glwe, ggsw)
	CHECK_CALL(vmp_apply_dft_p(module, result_dft, ncols, glwe->vec, nrows, N, ggsw_pmat, nrows, ncols),
		"vmp_apply_dft_p failed in ggsw_external_product");

	// Computes the bivGGSW ciphertext out of the DFT domain
	CHECK_CALL(vec_znx_idft_p(module, result->vec, ncols, result_dft, ncols), 
		"vec_znx_idft_p failed in ggsw_external_product");

	status = 0;

cleanup:
	free(result_dft);
	free(ggsw_pmat);

	return 0;
}

//! bivGGSW IN DFT PART (begin)

int glwe_secret_demasking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params_glwe, PolyBiv* result,
                                       const GGSWSecretKeyDFT* sk_dft, const VecBivDFT* glwe_vec_dft)
{
	int status = -1;

	// Point to the bivGLWE ciphertext out of the DFT domain
	VecBiv* glwe_vec = malloc(glwe_bytes(params_glwe));
	CHECK_ALLOC(glwe_vec, "glwe_vec's malloc failed in glwe_secret_demasking_ggsw_lib_dft");

	// Computes the bivGLWE ciphertext out of the DFT domain
	CHECK_CALL(vec_znx_idft_p(module, glwe_vec, glwe_size(params_glwe), glwe_vec_dft, glwe_size(params_glwe)),
		"vec_znx_idft_p failed in glwe_secret_demasking_ggsw_lib_dft");

	// Computes the phase in Zn[X,Y]
	CHECK_CALL(glwe_secret_demasking_ggsw_lib(module, params_glwe, result, sk_dft, glwe_vec), 
		"glwe_secret_demasking_ggsw_lib failed glwe_secret_demasking_ggsw_lib_dft");

	status = 0;

cleanup:
	free(glwe_vec);
	
	return 0;
}

int glwe_secret_masking_ggsw_lib_dft(const MODULE* module, const GLWECtParams* params_glwe, VecBivDFT* result_dft,
                                     const GGSWSecretKeyDFT* sk_dft, const PolyBivDFT* phase_dft)
{
	int status = -1;

	// bivGLWE parameters
	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t kappa   = params_glwe->kappa;
	uint64_t n_limbs = params_glwe->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// Variables 
	VecBiv* glwe_vec = NULL;
	PolyBiv* acc = NULL;
	PolyBivDFT* as_j_dft = NULL;
	PolyBiv* as_j = NULL;

	// Temporary bivGLWE ciphertext
	glwe_vec = malloc(N * l * (k + 1) * sizeof(int64_t));
	CHECK_ALLOC(glwe_vec, "glwe_vec's malloc failed");

	// Draws uniformly the ajs'
	CHECK_CALL(uniform_random_vec(k * N, glwe_vec, l, (k + 1) * N, kappa), 
		"uniform_random_vec failed in glwe_secret_masking_ggsw_lib_dft");

	// acc_(j+1) = acc_j + (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
	acc = calloc(N * l, sizeof(double));
	CHECK_ALLOC(acc, "acc's calloc failed in glwe_secret_masking_ggsw_lib_dft");

	// Point to resVec_j_dft = (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
	as_j_dft = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(as_j_dft, "as_j_dft's malloc failed in glwe_secret_masking_ggsw_lib_dft");

	// Point to resVec_j in Zn[X,Y] space
	as_j = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(as_j, "as_j's malloc failed in glwe_secret_masking_ggsw_lib_dft");

	// Computes Sum_j{0,k-1}[resVec_j]
	for (uint64_t j = 0; j < k; j++)
	{
		// The j-ème component of the secret key sk_dft
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];

		// Computes resVec_j_dft = (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, glwe_vec + j * N, l, (k + 1) * N);

		// Computes resVec_j in Zn[X,Y] space
		CHECK_CALL(vec_znx_idft_p(module, as_j, l, as_j_dft, l), "vec_znx_idft_p failed in glwe_secret_masking_ggsw_lib_dft");

		// And adds it to acc_j : acc_(j+1) = acc_j + resVec_j
		for (uint64_t p = 0; p < N * l; p++)
		{
			acc[p] += as_j[p];
		}
	}
	

	// The pointer to limb_0(b) in Zn[X,Y]
	PolyUniv* b_0_univ = glwe_vec + k * N;

	// For each i in {0,l} limb_i(b) = acc_i = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	// Then b is normalized
	CHECK_CALL(vec_znx_normalize_base2k_p(module, kappa, b_0_univ, l, (k + 1) * N, acc, l, N),
		"vec_normalize_base2k_p failed in glwe_secret_masking_ggsw_lib_dft");

	// Computes glwe_vec in the DFT domain
	vec_znx_dft_p(module, result_dft, l * (k + 1), glwe_vec, l * (k + 1), N);

	// Adds the phase (message with error) to bivGLWE(0), the result is a bivGLWE(m + err) ciphertext 
	for (uint64_t i = 1; i <= l; i++)
	{
		for (uint64_t p = 0; p < N; p++)
		{
			// Adds DFT(limb_i(phase)) to DFT(limb_i(b))
			result_dft[(i - 1) * (k + 1) * N + k * N + p] += phase_dft[(i - 1) * N + p];
		}
	}

cleanup: 
	free(as_j);
	free(as_j_dft);
	free(acc);
	free(glwe_vec);

	return 0;
}

int ggsw_secret_encrypt_dft(const MODULE* module, const GGSWCtParams* params_ggsw, GGSWCiphertextDFT* result_dft,
                            const GGSWSecretKeyDFT* sk_dft, const PolyUniv* m_univ)
{
	if (params_ggsw->k_tilde > params_ggsw->params_glwe->k)
		return log_perror("k_tilde should not be greater than k in ggsw_secret_encrypt_dft");

	int status = -1;

	const GLWECtParams* params_glwe = params_ggsw->params_glwe;

	// bivGLWE parameters
	uint64_t N       = params_glwe->N;
	uint64_t k       = params_glwe->k;
	uint64_t k_tilde = params_ggsw->k_tilde;

	// Variables
	PolyUnivDFT* m_univ_dft = NULL;   // DFT(msg)
	PolyUnivDFT* m_skj_univ_dft = NULL; // DFT(msg * sk_j)
	PolyUniv* m_skj_univ = NULL;        // -m*sk_j

	double* msk_univ_RnX = NULL;
	PolyBiv* glwe_biv_msg = NULL; // Temp space for -m * sk_j / 2^{...}
	PolyBivDFT* glwe_biv_msg_dft = NULL;


	// Point to DFT(msg)
	m_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_univ_dft, "m_univ_dft's malloc failed");
	
	// Computes DFT(m)
	vec_znx_dft_p(module, m_univ_dft, 1, m_univ, 1, N);

	// Point to DFT(m * sk_j)
	m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ_dft, "m_skj_univ_dft's malloc failed");

	// Point to -m * sk_j
	m_skj_univ = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(m_skj_univ, "m_skj_univ's alloc failed");

	// Point to (in univariate space) : -m * sk_j / 2^{kappa_tilde*i}, if j < k
	//                                            m / 2^{kappa_tilde*i}, if j = k
	msk_univ_RnX = malloc(poly_univ_bytes(params_glwe));
	CHECK_ALLOC(msk_univ_RnX, "m_skj_univ's alloc failed");

	// Point to (in bivariate space) : Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
	//                                           Dec_Kappa(m / 2^{kappa_tilde*i}) + err, if j = k
	glwe_biv_msg = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg, "phase's malloc failed");

	// We'll store DFT(-m * sk_j / 2^{kappa_tilde*i})
	glwe_biv_msg_dft = malloc(poly_biv_bytes(params_glwe));
	CHECK_ALLOC(glwe_biv_msg_dft, "phase_dft's malloc failed");

	for (uint64_t i = 1; i <= nb_partials(params_ggsw); i++)
	{
		for (uint64_t j = 0; j < k_tilde + 1; j++)
		{
			// Computes: Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
			//           Dec_Kappa(m / 2^{kappa_tilde*i}) + err,         if j = k
			// The precision of the decomposition is l

      if (j < k)
      {
        // Computes DFT(msg * sk_j)
        mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, m_univ_dft, 1);

        // Computes -DFT(msg * sk_j)
        // TODO: study if accelerable using AVX, do we need it in spqlios?
        for (uint64_t p = 0; p < N; p++) m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];

        // Inverse DFT to retreive -msg * sk_j
        CHECK_CALL(vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1), "vec_znx_idft_p failed in compute_phase_ij_dft");
      }
      for (uint64_t p = 0; p < N; p++)
        msk_univ_RnX[p] = ldexp( (j == k) ? (double)m_univ[p] : m_skj_univ[p], -params_ggsw->kappa_tilde * i);

      // Convert the result to a bivariate (base-2k) polynomial
      CHECK_CALL(univ_to_biv(params_glwe, glwe_biv_msg, msk_univ_RnX), "univ_to_biv failed in compute_phase_ij_dft");

      // Adds the error to get Dec_Kappa(m / 2^{kappa_tilde*}) + err
      // or                    Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err
      CHECK_CALL(add_bivariate_error(module, params_glwe, glwe_biv_msg, glwe_biv_msg), "add_error failed in compute_phase_ij_dft");

      // Permorm DFT to get the result in the DFT domain
      vec_znx_dft_p(module, glwe_biv_msg_dft, poly_biv_size(params_glwe), glwe_biv_msg, poly_biv_size(params_glwe), N);

      // Result destination
      VecBivDFT* glwe_vec_dft = ggsw_Sj_Yti_dft(params_ggsw, result_dft->mat, j, i);

      // Finally, mask/encrypt the above result to get a bivGLWE
      CHECK_CALL(glwe_secret_masking_ggsw_lib_dft(module, params_glwe, glwe_vec_dft, sk_dft, glwe_biv_msg_dft),
            "glwe_secret_masking_ggsw_lib_dft failed in ggsw_secret_encrypt_dft");
        }
      }

	status = 0;

cleanup:
	free(glwe_biv_msg_dft);
	free(glwe_biv_msg);
	free(msk_univ_RnX);
	free(m_skj_univ);
	free(m_skj_univ_dft);
	free(m_univ_dft);

	return 0;
}

int ggsw_external_product_dft(const MODULE* module,
                              GLWECiphertextDFT* result_dft,      // result
                              const GLWECiphertextDFT* glwe_dft,  // bivGLWE ciphertext
                              const GGSWCiphertextDFT* ggsw_dft   // bivGGSW ciphertext
)
{
	int status = -1;

	// Degree of chosen cyclotomic polynomial
	uint64_t N = result_dft->params->N;

	// The bivGLWE ciphertext glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// The bivGGSW ciphertext ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// As the result of the vector-matrix product glwe * ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ggsw_dft->params->n_limbs_tilde;
	uint64_t ncols = ggsw_dft->params->params_glwe->n_limbs;

	// Variables
	MatBiv* ggsw_mat = NULL;
	MatBivDFT* ggsw_pmat = NULL;

	// Point to the bivGGSW ciphertext out of the DFT domain
	ggsw_mat = malloc(ggsw_bytes(ggsw_dft->params));
	CHECK_ALLOC(ggsw_mat, "mat's malloc failed in ggsw_external_product_dft");
	
	// Computes the bivGGSW ciphertext out of the DFT domain
	CHECK_CALL(vec_znx_idft_p(module, ggsw_mat, nrows * ncols, ggsw_dft->mat, nrows * ncols),"vec_znx_idft_p failed in ggsw_external_product_dft");

	// Point to the bivGGSW ciphertext in the DFT domain
	ggsw_pmat = malloc(ggsw_bytes(ggsw_dft->params));
	CHECK_ALLOC(ggsw_pmat, "pmat's malloc failed in ggsw_external_product_dft");

	// Prepares the bivGGSW ciphertext in the DFT domain
	CHECK_CALL(vmp_prepare_contiguous_p(module, ggsw_pmat, ggsw_mat, nrows, ncols), "vmp_prepare_contiguous_p failed in ggsw_external_product_dft");

	// Computes ExternalProduct(glwe, ggsw)
	CHECK_CALL(vmp_apply_dft_to_dft_p(module, result_dft->vec, ncols, glwe_dft->vec, nrows, ggsw_pmat, nrows, ncols), 
			  "vmp_apply_dft_to_dft_p failed in ggsw_external_product_dft");


	status = 0;

cleanup:
	free(ggsw_pmat);
	free(ggsw_mat);

	return status;
}
