# bivGLWE Documentation

This file documents the internal representation of a bivariate bivGLWE ciphertext used in the library. It focuses on the mathematical structure and memory layout rather than the API.

## Notations
### Polynomials

Let P(X) be a polynomial in Z_n[X]. Another representation is (P_0, ... $cdot$)   


/**
 * @page glwe_bivariate_ciphertext Bivariate bivGLWE Ciphertext
 *
 * @section glwe_bivariate_overview glwe_bivariate_overview
 *
 * A bivGLWE ciphertext \f$\in (\mathbb{Z}_n[X])^k\f$ encodes a polynomial in @p X 
 * with coefficients in \f$\mathbb{Z}\f$.
 *
 * A bivariate bivGLWE ciphertext \f$\in (\mathbb{Z}_n[X])^k\f$ encodes a polynomial 
 * in @p X with coefficients in \f$\mathbb{Z}\f$.
 *
 *
 * It generalizes the classical bivGLWE ciphertext by organizing the
 * polynomial coefficients along a second dimension @p Y, allowing
 * structured decompositions and gadget-based constructions.
 *
 * The ciphertext is represented as a flattened array of coefficients
 * corresponding to the polynomials:
 *
 * @verbatim
 * (a_0(Y), a_1(Y), ..., a_{k-1}(Y), b(Y))
 * @endverbatim
 *
 * where each @p a_i(Y) and @p b(Y) are univariate polynomials in @p Y.
 *
 * @section glwe_bivariate_structure Internal structure
 *
 * The ciphertext is stored in a single contiguous array @p ct of
 * length @p n_limbs.
 *
 * @verbatim
 * ct = [ a_0[0], a_1[0], ..., a_k[0] ]
 *        a_0[1], a_1[1], ..., a_k[1],
 *        ...
 *      
 * @endverbatim
 *
 * The index mapping is defined as follows:
 *
 * - Let @p k be the bivGLWE dimension
 * - Let @p n_limbs be the total number of limbs
 * - Let:
 *   - @p l_a = floor((n_limbs + 1) / (k + 1))
 *   - @p l_b = n_limbs - k * l_a
 *
 * For a given position @p p in @p ct:
 *
 * @verbatim
 * i = p % (k + 1)          // polynomial index
 * j = floor(p / (k + 1))  // Y-degree index
 * @endverbatim
 *
 * The coefficient corresponds to:
 *
 * @verbatim
 * ct[p] = a_i[j]
 * @endverbatim
 *
 * where:
 * - @p i in [0, k-1] refers to polynomial @p a_i
 * - @p i = k refers to the polynomial @p b
 *
 * @section glwe_bivariate_polynomials Polynomial interpretation
 *
 * The ciphertext encodes the following polynomials in @p Y:
 *
 * @verbatim
 * a_0(Y) = a_0[0] + a_0[1] Y + ... + a_0[l_a - 1] Y^{l_a - 1}
 * a_1(Y) = a_1[0] + a_1[1] Y + ... + a_1[l_a - 1] Y^{l_a - 1}
 * ...
 * a_{k-1}(Y) = a_{k-1}[0] + ... + a_{k-1}[l_a - 1] Y^{l_a - 1}
 *
 * b(Y) = b[0] + b[1] Y + ... + b[l_b - 1] Y^{l_b - 1}
 * @endverbatim
 *
 * The polynomial @p b(Y) may have a different degree bound than the
 * @p a_i(Y) polynomials.
 *
 * @section glwe_bivariate_precision Precision and decomposition
 *
 * The ciphertext precision is controlled by the parameter:
 *
 * @verbatim
 * L = l_b * kappa
 * @endverbatim
 *
 * where @p kappa is the base decomposition parameter.
 *
 * This structure is particularly suited for gadget decompositions
 * and bivariate constructions appearing in advanced homomorphic
 * encryption schemes.
 *
 * @section glwe_bivariate_summary Summary
 *
 * A bivariate bivGLWE ciphertext can be viewed as:
 *
 * - A bivGLWE ciphertext whose coefficients are polynomials in @p Y
 * - A matrix-like organization flattened into a 1D array
 * - A flexible structure enabling precision control and gadget-based
 *   operations
 */
