#include "ggsw.h"

#include <string.h>

#include "glwe_ciphertext.h"
#include "logger.h"
#include "math.h"
#include "rng.h"
#include "spqlios_alias.h"

//! GGSW PART (begin)

void add_error(MODULE* module, GLWECtParams* params, PolyBiv* res, PolyBiv* phase)
{
	// Compute a random error in DFT space
	PolyBiv* err = new_normal_random_biv_poly(module, params);

	// Add the error in DFT space
	add_biv_poly(params, res, params->N, phase, params->N, err, params->N);

	free(err);
}

int glwe_secret_demasking_ggsw_lib(MODULE* module, GLWECtParams* params, PolyBiv* phase, GGSWSecretKeyDFT* sk_dft,
                                   VecBiv* ct)
{
	// GLWE parameters
	uint64_t N   = params->N;
	uint64_t k   = params->k;
	uint64_t l   = poly_biv_size(params);

	PolyBiv* acc = calloc(N * l, sizeof(int64_t));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_demasking_ggsw_lib.") < 0) return -1;

	// Will point to DFT(sk_j * a_j)
	PolyBivDFT* as_j_dft = calloc(2 * poly_biv_coef_number_dft(params), sizeof(double));
	if (log_is_null(as_j_dft, "as_j_dft's calloc failed in glwe_secret_demasking_ggsw_lib.") < 0) {
		free(acc);
		return -1;
	}

	// Will point to sk_j * a_j
	PolyBiv* as_j = calloc(poly_biv_coef_number(params), sizeof(int64_t));
	if (log_is_null(as_j, "as_j's calloc failed in glwe_secret_demasking_ggsw_lib.") < 0) {
		free(acc);
		free(as_j_dft);
		return -1;
	}

	// Computes acc = -Sum_j{0,k-1}[sk_j * a_j]
	for (int64_t j = 0; j < k; j++) {
		// The j-ème component of resp. the secret key and the bivGLWE ciphertext
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];
		PolyUniv* a_j              = ct + j * N;

		// Computes DFT(sk_j * a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, a_j, l, (k + 1) * N);

		// Computes sk_j * a_j
		vec_znx_idft_p(module, as_j, l, as_j_dft, l);

		// And subs it to acc
		for (int64_t p = 0; p < N * l; p++) {
			acc[p] -= as_j[p];
		}
	}
	free(as_j_dft);
	free(as_j);

	// Computes acc = b - Sum_j{0,k-1}[sk_j * a_j]
	int64_t* b = ct + k * N;
	add_biv_poly(params, acc, N, b, (k + 1) * N, acc, N);

	vec_znx_normalize_base2k_p(module, params->kappa, phase, l, N, acc, l, N);

	free(acc);

	return 0;
}

int glwe_secret_masking_ggsw_lib(MODULE* module, GLWECtParams* params, VecBiv* res_ct, GGSWSecretKeyDFT* sk_dft,
                                 PolyBiv* phase)
{
	uint64_t N       = params->N;
	uint64_t k       = params->k;
	uint64_t kappa   = params->kappa;
	uint64_t n_limbs = params->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	if (inplace_uniform_random_vec(k * N, res_ct, l, (k + 1) * N, kappa) < 0) return -1;

	// acc_(j+1) = acc_j + (sk_j * limb_1(a_j) , ... , sk_j * limb_l(a_j))
	PolyBiv* acc = calloc(N * l, sizeof(double));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_masking_ggsw_lib.") < 0) return -1;

	// Will point to DFT(sk_j) * DFT(a_j)
	PolyBivDFT* as_j_dft = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j_dft, "as_j_dft's calloc failed in glwe_secret_masking_ggsw_lib.") < 0) {
		free(acc);
		return -1;
	}

	// Will point to sk_j * a_j
	PolyBiv* as_j = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j, "as_j's calloc failed in glwe_secret_masking_ggsw_lib") < 0) {
		free(acc);
		free(as_j_dft);
		return -1;
	}

	// Computes Sum_j{0,k-1}[sk_j * a_j]
	for (int64_t j = 0; j < k; j++) {
		// The j-ème component of the secret key sk_dft
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];

		// Computes DFT(sk_j) * DFT(a_j)
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, res_ct + j * N, l, (k + 1) * N);

		// Computes sk_j * a_j
		vec_znx_idft_p(module, as_j, l, as_j_dft, l);

		// And adds it to acc_j
		for (int64_t p = 0; p < N * l; p++) {
			acc[p] += as_j[p];
		}
	}
	free(as_j_dft);
	free(as_j);

	// Add the phase to acc
	for (int64_t i = 1; i <= l; i++) {
		for (int64_t p = 0; p < N; p++) {
			acc[(i - 1) * N + p] += phase[(i - 1) * N + p];
		}
	}

	// The pointer to limb_0(b)
	PolyBiv* b_0 = res_ct + k * N;

	// For each i in {0,l} limb_i(b) = limb_i(acc) = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	vec_znx_normalize_base2k_p(module, kappa, b_0, l, (k + 1) * N, acc, l, N);

	free(acc);

	return 0;
}

int compute_phase_ij(MODULE* module, GGSWCtParams* params_ggsw, GGSWSecretKeyDFT* sk_dft, PolyUniv* msg_univ,
                     PolyUnivDFT* msg_univ_dft, PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft, PolyBiv* phase,
                     PolyUnivRnX* phase_univ_RnX, int64_t i, int64_t j)
{
	// GLWE parameters
	GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                = params_glwe->N;

	if (j < params_glwe->k) {
		// Computes DFT(msg * sk_j)
		mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);

		// Computes -DFT(msg * sk_j)
		for (int64_t p = 0; p < N; p++) {
			m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
		}

		// Computes -msg * sk_j
		vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1);

		// Computes -msg * sk_j / 2^{kappa_tilde*i}
		for (int64_t p = 0; p < N; p++) {
			phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -params_ggsw->kappa_tilde * i);
		}

		// Compute the base-2^kappa decomposition of : -m * sk_j / 2^{kappa_tilde*i}
		if (univ_to_biv(params_glwe, phase, phase_univ_RnX) < 0)
			return log_perror("univ_to_biv failed in compute_phase_ij.");

		// Computes the phase Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err
		add_error(module, params_glwe, phase, phase);
	} else {
		// Computes msg / 2^{kappa_tilde*i}
		for (int64_t p = 0; p < N; p++) {
			phase_univ_RnX[p] = ldexp((double)msg_univ[p], -params_ggsw->kappa_tilde * i);
		}

		// Compute the base-2^kappa decomposition of : m / 2^{kappa_tilde*i}
		if (univ_to_biv(params_glwe, phase, phase_univ_RnX) < 0)
			return log_perror("univ_to_biv failed in compute_phase_ij.");

		// Computes the phase Dec_Kappa(m / 2^{kappa_tilde*i}) + err
		add_error(module, params_glwe, phase, phase);
	}
}

int ggsw_secret_encrypt(MODULE* module, GGSWCtParams* params_ggsw, GGSWCiphertext* ct_ggsw, GGSWSecretKeyDFT* sk_dft,
                        PolyUniv* msg_univ)
{
	if (params_ggsw->k_tilde > params_ggsw->params_glwe->k)
		return log_perror("k_tilde should not be greater than k in ggsw_secret_encrypt.");

	// GLWE parameters
	GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                = params_glwe->N;
	uint64_t k                = params_glwe->k;
	uint64_t k_tilde          = params_ggsw->k_tilde;

	// Computes DFT(msg)
	PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(msg_univ_dft, "msg_univ_dft's malloc failed in ggsw_secret_encrypt.") < 0) return -1;
	vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, N);

	// Will store DFT(msg * sk_j)
	PolyUnivDFT* m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(m_skj_univ_dft, "m_skj_univ_dft's malloc failed in ggsw_secret_encrypt.") < 0) {
		free(msg_univ_dft);
		return -1;
	}

	// Will store -msg * sk_j
	PolyUniv* m_skj_univ = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(m_skj_univ, "m_skj_univ's malloc failed in ggsw_secret_encrypt.") < 0) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		return -1;
	}

	// Will store (in univariate space) : -m * sk_j / 2^{kappa_tilde*i}, if j < k
	//                                            m / 2^{kappa_tilde*i}, if j = k
	double* phase_univ_RnX = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(phase_univ_RnX, "phase_univ_RnX's malloc failed in ggsw_secret_encrypt.") < 0) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		free(m_skj_univ);
		return -1;
	}

	// Will store (in bivariate space) : Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
	//                                           Dec_Kappa(m / 2^{kappa_tilde*i}) + err, if j = k
	PolyBiv* phase = malloc(poly_biv_bytes(params_glwe));
	if (log_is_null(phase, "phase's malloc failed in ggsw_secret_encrypt.")) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		free(m_skj_univ);
		free(phase_univ_RnX);
		return -1;
	}

	for (int64_t i = 1; i <= nb_partials(params_ggsw); i++) {
		for (int64_t j = 0; j < k_tilde + 1; j++) {
			// Computes the the bivariate phase : Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
			//                                            Dec_Kappa(m / 2^{kappa_tilde*i}) + err, if j = k
			// The precision of the decomposition is l
			if (compute_phase_ij(module, params_ggsw, sk_dft, msg_univ, msg_univ_dft, m_skj_univ, m_skj_univ_dft, phase,
			                     phase_univ_RnX, i, j) < 0) {
				free(phase);
				free(msg_univ_dft);
				free(m_skj_univ_dft);
				free(m_skj_univ);
				free(phase_univ_RnX);
				return log_perror("compute_phase_ij failed in ggsw_secret_encrypt.");
				;
			}

			// The pointer to : bivGLWE(-m * sk_j / 2^{kappa_tilde*i}), if j < k
			//                         bivGLWE( m / 2^{kappa_tilde*i}), if j = k
			VecBiv* ct_glwe = ggsw_Sj_Yti(params_ggsw, ct_ggsw->mat, j, i);

			// Computes : bivGLWE(-m * sk_j / 2^{kappa_tilde*i}), if j < k
			//                   bivGLWE( m / 2^{kappa_tilde*i}), if j = k
			if (glwe_secret_masking_ggsw_lib(module, params_glwe, ct_glwe, sk_dft, phase) < 0) {
				free(phase);
				free(msg_univ_dft);
				free(m_skj_univ_dft);
				free(m_skj_univ);
				free(phase_univ_RnX);
				return log_perror("glwe_secret_masking_ggsw_lib failed in ggsw_secret_encrypt.");
				;
			}
		}
	}

	free(phase);
	free(msg_univ_dft);
	free(m_skj_univ_dft);
	free(m_skj_univ);
	free(phase_univ_RnX);

	return 0;
}

int ggsw_external_product(MODULE* module,
                          GLWECiphertext* res,      // result
                          GLWECiphertext* ct_glwe,  // GLWE ciphertext
                          GGSWCiphertext* ct_ggsw   // GGSW ciphertext
)
{
	uint64_t N = res->params->N;

	// The bivGGSW ciphertext ct_ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// The bivGLWE ciphertext ct_glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// As the result of the vector-matrix product ct_glwe * ct_ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ct_ggsw->params->n_limbs_tilde;
	uint64_t ncols = ct_ggsw->params->params_glwe->n_limbs;

	// Prepares the GGSW ciphertext in DFT space
	MatBivDFT* mat_dft = malloc(ggsw_bytes(ct_ggsw->params));
	if (log_is_null(mat_dft, "mat_dft's malloc failed in ggsw_external_product.") < 0) return -1;

	vmp_prepare_contiguous_p(module, mat_dft, ct_ggsw->mat, nrows, ncols);

	// The pointer to ExternalProduct(ct_glwe, ct_ggsw)
	VecBivDFT* result = malloc(glwe_bytes(ct_ggsw->params->params_glwe));
	if (log_is_null(result, "result's malloc failed in ggsw_external_product.") < 0) {
		free(mat_dft);
		return -1;
	}

	// Computes ExternalProduct(ct_glwe, ct_ggsw)
	vmp_apply_dft_p(module, result, ncols, ct_glwe->vec, nrows, N, mat_dft, nrows, ncols);
	vec_znx_idft_p(module, res->vec, ncols, result, ncols);

	free(mat_dft);
	free(result);
}

//! GGSW IN DFT PART (begin)

void add_error_dft(MODULE* module, GLWECtParams* params, PolyBivDFT* res_dft, PolyBivDFT* phase_dft)
{
	// Compute a random error in DFT space
	PolyBivDFT* err_dft = new_normal_random_biv_poly_dft(module, params);

	// Add the error in DFT space
	add_biv_poly_dft(params, phase_dft, params->N, phase_dft, params->N, err_dft, params->N);
}

int glwe_secret_demasking_ggsw_lib_dft(MODULE* module, GLWECtParams* params, PolyBiv* phase, GGSWSecretKeyDFT* sk_dft,
                                       VecBivDFT* ct_dft)
{
	// Computes the ciphertext out of DFT space
	VecBiv* ct = malloc(glwe_bytes(params));
	if (log_is_null(ct, "ct's malloc failed in glwe_secret_demasking_ggsw_lib_dft.") < 0) return -1;
	vec_znx_idft_p(module, ct, glwe_size(params), ct_dft, glwe_size(params));

	if (glwe_secret_demasking_ggsw_lib(module, params, phase, sk_dft, ct) < 0) {
		free(ct);
		return log_perror("glwe_secret_demasking_ggsw_lib failed glwe_secret_demasking_ggsw_lib_dft.");
	}

	free(ct);
	return 0;
}

int glwe_secret_masking_ggsw_lib_dft(MODULE* module, GLWECtParams* params, VecBivDFT* res_dft, GGSWSecretKeyDFT* sk_dft,
                                     PolyBivDFT* phase_dft)
{
	uint64_t N       = params->N;
	uint64_t k       = params->k;
	uint64_t kappa   = params->kappa;
	uint64_t n_limbs = params->n_limbs;
	uint64_t l       = n_limbs / (k + 1);

	// Temporary bivGLWE ciphertext
	VecBiv* tmp_ct = malloc(N * l * (k + 1) * sizeof(int64_t));
	if (log_is_null(tmp_ct, "tmp_ct's malloc failed.") < 0) return -1;

	if (inplace_uniform_random_vec(k * N, tmp_ct, l, (k + 1) * N, kappa) < 0) {
		free(tmp_ct);
		return log_perror("inplace_uniform_random_vec failed in glwe_secret_masking_ggsw_lib_dft.");
	}

	// acc_(j+1) = acc_j + (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
	PolyBiv* acc = calloc(N * l, sizeof(double));
	if (log_is_null(acc, "acc's calloc failed in glwe_secret_masking_ggsw_lib_dft.") < 0) return -1;

	// Will point to resVec_j_dft = (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
	PolyBivDFT* as_j_dft = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j_dft, "as_j_dft's malloc failed in glwe_secret_masking_ggsw_lib_dft.") < 0) {
		free(acc);
		return -1;
	}

	// Will point to resVec_j in Zn[X,Y] space
	PolyBiv* as_j = malloc(poly_biv_bytes(params));
	if (log_is_null(as_j, "as_j's malloc failed in glwe_secret_masking_ggsw_lib_dft") < 0) {
		free(acc);
		free(as_j_dft);
		return -1;
	}

	// Computes Sum_j{0,k-1}[resVec_j]
	for (int64_t j = 0; j < k; j++) {
		// The j-ème component of the secret key sk_dft
		PolyUnivDFT* sk_j_univ_dft = sk_dft->values[j];

		// Computes resVec_j_dft = (DFT(sk_j) * limb_1(a_j) , ... , DFT(sk_j) * limb_l(a_j))
		svp_apply_dft_p(module, as_j_dft, l, sk_j_univ_dft, tmp_ct + j * N, l, (k + 1) * N);

		// Computes resVec_j in Zn[X,Y] space
		vec_znx_idft_p(module, as_j, l, as_j_dft, l);

		// And adds it to acc_j : acc_(j+1) = acc_j + resVec_j
		for (int64_t p = 0; p < N * l; p++) {
			acc[p] += as_j[p];
		}
	}
	free(as_j);
	free(as_j_dft);

	// The pointer to limb_0(b) in Zn[X,Y]
	PolyUniv* b_0_univ = tmp_ct + k * N;

	// For each i in {0,l} limb_i(b) = acc_i = Sum_j{0,k-1}[sk_j * limb_i(a_j)]
	// Then b is normalized
	vec_znx_normalize_base2k_p(module, kappa, b_0_univ, l, (k + 1) * N, acc, l, N);

	// Computes tmp_ct in DFT space
	vec_znx_dft_p(module, res_dft, l * (k + 1), tmp_ct, l * (k + 1), N);

	// Adds the phase (message with error) to bivGLWE(0), the result is a ct of bivGLWE(m + e)
	for (int64_t i = 1; i <= l; i++) {
		for (int64_t p = 0; p < N; p++) {
			// Adds DFT(limb_i(phase)) to DFT(limb_i(b))
			res_dft[(i - 1) * (k + 1) * N + k * N + p] =
			    res_dft[(i - 1) * (k + 1) * N + k * N + p] + phase_dft[(i - 1) * N + p];
		}
	}

	free(acc);
	free(tmp_ct);

	return 0;
}

int compute_phase_ij_dft(MODULE* module, GGSWCtParams* params_ggsw, GGSWSecretKeyDFT* sk_dft, PolyUniv* msg_univ,
                         PolyUnivDFT* msg_univ_dft, PolyUniv* m_skj_univ, PolyUnivDFT* m_skj_univ_dft,
                         PolyBivDFT* phase_dft, PolyBiv* phase, PolyUnivRnX* phase_univ_RnX, int64_t i, int64_t j)
{
	// GLWE parameters
	GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                = params_glwe->N;

	if (j < params_glwe->k) {
		// Computes DFT(msg * sk_j)
		mult_vec_znx_dft(module, m_skj_univ_dft, 1, sk_dft->values[j], 1, msg_univ_dft, 1);

		// Computes -DFT(msg * sk_j)
		for (int64_t p = 0; p < N; p++) {
			m_skj_univ_dft[p] = -1 * m_skj_univ_dft[p];
		}

		// Computes -msg * sk_j
		vec_znx_idft_p(module, m_skj_univ, 1, m_skj_univ_dft, 1);

		// Computes -msg * sk_j / 2^{kappa_tilde*i}
		for (int64_t p = 0; p < N; p++) {
			phase_univ_RnX[p] = ldexp((double)m_skj_univ[p], -params_ggsw->kappa_tilde * i);
		}

		// Compute the base-2^kappa decomposition of : -m * sk_j / 2^{kappa_tilde*i}
		if (univ_to_biv(params_glwe, phase, phase_univ_RnX) < 0)
			return log_perror("univ_to_biv failed in compute_phase_ij_dft.");

		// Computes the phase Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err
		add_error(module, params_glwe, phase, phase);

		vec_znx_dft_p(module, phase_dft, poly_biv_size(params_glwe), phase, poly_biv_size(params_glwe), N);

		return 0;
	} else {
		// Computes msg / 2^{kappa_tilde*i}
		for (int64_t p = 0; p < N; p++) {
			phase_univ_RnX[p] = ldexp((double)msg_univ[p], -params_ggsw->kappa_tilde * i);
		}

		// Compute the base-2^kappa decomposition of : m / 2^{kappa_tilde*i}
		if (univ_to_biv(params_glwe, phase, phase_univ_RnX) < 0)
			return log_perror("univ_to_biv failed in compute_phase_ij_dft.");

		// Computes the phase Dec_Kappa(m / 2^{kappa_tilde*}) + err
		add_error(module, params_glwe, phase, phase);

		vec_znx_dft_p(module, phase_dft, poly_biv_size(params_glwe), phase, poly_biv_size(params_glwe), N);

		return 0;
	}
}

int ggsw_secret_encrypt_dft(MODULE* module, GGSWCtParams* params_ggsw, GGSWCiphertextDFT* ct_ggsw_dft,
                            GGSWSecretKeyDFT* sk_dft, PolyUniv* msg_univ)
{
	if (params_ggsw->k_tilde > params_ggsw->params_glwe->k)
		return log_perror("k_tilde should not be greater than k in ggsw_secret_encrypt_dft.");

	// GLWE parameters
	GLWECtParams* params_glwe = params_ggsw->params_glwe;
	uint64_t N                = params_glwe->N;
	uint64_t k                = params_glwe->k;
	uint64_t k_tilde          = params_ggsw->k_tilde;

	// Computes DFT(msg)
	PolyUnivDFT* msg_univ_dft = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(msg_univ_dft, "msg_univ_dft's malloc failed.") < 0) return -1;
	vec_znx_dft_p(module, msg_univ_dft, 1, msg_univ, 1, N);

	// Will store DFT(msg * sk_j)
	PolyUnivDFT* m_skj_univ_dft = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(m_skj_univ_dft, "m_skj_univ_dft's malloc failed.") < 0) {
		free(msg_univ_dft);
		return -1;
	}

	// Will store -msg * sk_j
	PolyUniv* m_skj_univ = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(m_skj_univ, "m_skj_univ's alloc failed.") < 0) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		return -1;
	}

	// Will store (in univariate space) : -m * sk_j / 2^{kappa_tilde*i}, if j < k
	//                                            m / 2^{kappa_tilde*i}, if j = k
	double* phase_univ_RnX = malloc(poly_univ_bytes(params_glwe));
	if (log_is_null(phase_univ_RnX, "m_skj_univ's alloc failed.") < 0) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		free(m_skj_univ);
		return -1;
	}

	// Will store (in bivariate space) : Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
	//                                           Dec_Kappa(m / 2^{kappa_tilde*i}) + err, if j = k
	PolyBiv* phase = malloc(poly_biv_bytes(params_glwe));
	if (log_is_null(phase, "phase's malloc failed.") < 0) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		free(m_skj_univ);
		free(phase_univ_RnX);
		return -1;
	}

	// We'll store DFT(-m * sk_j / 2^{kappa_tilde*i})
	PolyBivDFT* phase_dft = malloc(poly_biv_bytes(params_glwe));
	if (log_is_null(phase_dft, "phase_dft's malloc failed.") < 0) {
		free(msg_univ_dft);
		free(m_skj_univ_dft);
		free(m_skj_univ);
		free(phase_univ_RnX);
		free(phase);
		return -1;
	}

	for (int64_t i = 1; i <= nb_partials(params_ggsw); i++) {
		for (int64_t j = 0; j < k_tilde + 1; j++) {
			// Computes the the bivariate phase : Dec_Kappa(-m * sk_j / 2^{kappa_tilde*i}) + err, if j < k
			//                                            Dec_Kappa(m / 2^{kappa_tilde*i}) + err, if j = k
			// The precision of the decomposition is l
			if (compute_phase_ij_dft(module, params_ggsw, sk_dft, msg_univ, msg_univ_dft, m_skj_univ, m_skj_univ_dft,
			                         phase_dft, phase, phase_univ_RnX, i, j) < 0) {
				free(phase);
				free(msg_univ_dft);
				free(m_skj_univ_dft);
				free(m_skj_univ);
				free(phase_univ_RnX);
				return log_perror("compute_phase_ij_dft failed in ggsw_secret_encrypt_dft.");
			}

			// The pointer in DFT space to : bivGLWE(-m * sk_j / 2^{kappa_tilde*i}), if j < k
			//                               bivGLWE( m / 2^{kappa_tilde*i}), if j = k
			VecBivDFT* ct_glwe = ggsw_Sj_Yti_dft(params_ggsw, ct_ggsw_dft->mat, j, i);

			// Computes in DFT space: bivGLWE(-m * sk_j / 2^{kappa_tilde*i}), if j < k
			//                        bivGLWE( m / 2^{kappa_tilde*i}), if j = k
			if (glwe_secret_masking_ggsw_lib_dft(module, params_glwe, ct_glwe, sk_dft, phase_dft) < 0) {
				free(phase);
				free(m_skj_univ);
				free(phase_univ_RnX);
				free(phase_dft);
				free(msg_univ_dft);
				free(m_skj_univ_dft);
				return log_perror("glwe_secret_demasking_ggsw_lib_dft failed in ggsw_secret_encrypt_dft.");
			}
		}
	}

	free(phase);
	free(m_skj_univ);
	free(phase_univ_RnX);
	free(phase_dft);
	free(m_skj_univ_dft);
	free(msg_univ_dft);

	return 0;
}

int ggsw_external_product_dft(MODULE* module,
                              GLWECiphertextDFT* res_dft,      // result
                              GLWECiphertextDFT* ct_glwe_dft,  // GLWE ciphertext
                              GGSWCiphertextDFT* ct_ggsw_dft   // GGSW ciphertext
)
{
	uint64_t N = res_dft->params->N;

	// The bivGLWE ciphertext ct_glwe is a prepared vector in Vec(Zn[X]) of size n_limbs_tilde
	// The bivGGSW ciphertext ct_ggsw is a prepared matrix in Mat(Zn[X]) of size n_limbs_tilde * n_limbs
	// As the result of the vector-matrix product ct_glwe * ct_ggsw,
	// the bivGLWE ciphertext res is a prepared vector in Vec(Zn[X]) of size n_limbs
	uint64_t nrows = ct_ggsw_dft->params->n_limbs_tilde;
	uint64_t ncols = ct_ggsw_dft->params->params_glwe->n_limbs;

	// Computes the GGSW ciphertext out of DFT space
	MatBiv* mat = malloc(ggsw_bytes(ct_ggsw_dft->params));
	if (log_is_null(mat, "mat's malloc failed in ggsw_external_product_dft.") < 0) return -1;
	vec_znx_idft_p(module, mat, nrows * ncols, ct_ggsw_dft->mat, nrows * ncols);

	// Prepares the GGSW ciphertext in DFT space
	MatBivDFT* pmat = malloc(ggsw_bytes(ct_ggsw_dft->params));
	if (log_is_null(pmat, "pmat's malloc failed in ggsw_external_product_dft.") < 0) {
		free(mat);
		return -1;
	}
	vmp_prepare_contiguous_p(module, pmat, mat, nrows, ncols);

	// Computes ExternalProduct(ct_glwe, ct_ggsw)
	vmp_apply_dft_to_dft_p(module, res_dft->pvec, ncols, ct_glwe_dft->pvec, nrows, pmat, nrows, ncols);

	free(mat);
	free(pmat);
}
