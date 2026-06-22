# Polynomials in prividema-fhe

This file explains the memory layout and considerations of the different representations of polynomials in Prividema-fhe.

One is expected to have read and be familiar with the [main paper used in this library](https://eprint.iacr.org/2023/771).
The mathematical concepts used will not be explained in this documentation. Rather, the implementation particularities will be.

In particular, \N, \K, \f$\ell\f$, etc. will be used here as in the paper.

## Univariate polynomial representations

All univariate representations are implemented by the means of a dynamic array, i.e., a virtually contiguous region of memory in the heap.

Let \f$ P(x) = \sum_{i = 0}^{N-1} p_i x^i \f$ be a polynomial.
If ``p_univ`` is the pointer (address) of the polynomial, element \f$ p_i \f$ will be found at
``p_univ + sizeof(*p_univ)*i`` (in bytes) or  `` p_univ + i `` (in C pointer arithmetic).

### \ZnX polynomials

\ZnX polynomials are represented using a signed 64bit integer (`int64_t`).
Overflow is undefined behaviour in C as well as in the library.

### \RnX polynomials

The \RnX representation might be the most intuitive one when elements in the torus could be involved.
Elements are represented using the type ``double``, which is assumed to correspond to
[IEEE-754 binary64](https://en.wikipedia.org/wiki/Double-precision_floating-point_format) encoding.
While strict compliance with the standard is not needed, the library might assume that
the memory encoding of `double` values is that of the standard.

As a summary `binary64` floating point numbers have 53 bits of significant precision, of which one bit is implicit,
and the rest of bits are dedicated to the sign and (binary) exponent. This gives a relative precision ([machine epsilon](https://en.wikipedia.org/wiki/Machine_epsilon))
of \f$ 2^{-53} \approx 1.1 \cdot 10^{-16} \f$, and a range of magnitudes from approximately \f$ 10 ^ {-308} \f$ to \f$ 10 ^ {308} \f$.

> [!WARNING]
> When doing an operation (eg., a sum) with two floating point numbers, one should assume that
> the result might have an absolute error of at least \f$ 2^{-53} \f$ times the largest number involved (inputs and outputs).
>
> This is very important if one plans on adding FHE noise in \RnX: all noise below approximately \f$ 1.1 \cdot 10^{-16}\f$ times
> the value to which it is being added (be it \f$ as\f$, \f$ m \f$ or \f$ as + m \f$) will be discarded!!

The \RnX representation can be used for elements in the torus as well as general real numbers.
Performance-wise, it is one of the least efficient representations due to being floating point.

Moreover, it is prone to underflow errors when used for torus elements (warning box above):
While in theory elements \f$ -2^{-54} \f$ and \f$ 1 -2^{-54} \f$ represent the same element of the torus,
due to floating point precision limitations, the first element will be encoded properly but the second will
be encoded as a 1 or equivalently in the torus, a 0!

### \TnX fixed-point representation

The fixed-point \TnX representation (called simply `tnx` in the codebase) can be used to represent elements in the torus
in a more efficient way.

Technically, it stores elements as the numerator of a fraction with implicit denominator \f$ 2^{-64} \f$.
Since the elements are stored as a 2-complement integers, torus arithmetic is simply its native 2-complement arithmetic.
The values can be interpreted as signed or unsigned as desired, again thanks to 2-complement arithmetic.
Thus, this format is the best for performance. It also has more bits of precision available than floating point.
In contrast, however, it has a constant denominator, which means that an absolute error of up to \f$ 2^{-64} \f$
should be assumed.

## Bivariate polynomials

This section details the notable architectural details of the bivariate polynomial representation, specifically addressing precision constraints, overflow management, and the flattened memory allocation layout in the following subsections.

### A note on precision and overflows

Internally, bivariate polynomials use `int64_t` for their elements (thus have 64bits of precision).
If using a FFT backend for computation, \K will typically be no larger than 19 (due to FFT precision limitations), which means that
\f$ \|P\|_\infty \leq 2^{18}\f$.

The `int64_t` precision gives us a bound for either \N the degree of our polynomials or \K the
number of bits per limb of a base-\f$ 2^K \f$ decomposition regardless of the backend used, since we must be able to
accurately represent polynomial products.
The backend will probably impose tighter bounds in most cases, yet the following should be true in all cases:

\f[
  2 \cdot (K-1) + \log_2(N) < 63
\f]

which results from the following inequality, itself the result of imposing that all coefficients of the
result are able to be represented in a `int64_t`:

\f[
\| A \cdot B \|_\infty \leq \max_{i}  \sum_{j + k \equiv i \pmod{N}} \left| a_{j} \cdot b_{k} \right| \leq 2^{63} - 1
\f]

### Representation in memory {#bivariate-notation}

As stated, bivariate polynomial coefficients are stored as `int64_t` values.
For a single bivariate polynomial, the coefficients are stored in a limb-major fashion.
In other words, if we have decomposition of \f$A \f$ which is \f$ A_{biv}(X,Y) = \sum_{i=0}^{N-1} \sum_{j=1}^{\ell} a_{i,j} X^i Y^j  \f$ such that
\f$ A_{biv}(X, 2^{-K}) = A(X) \f$, we store coefficients as follows (flattened in a row-major order):

\f[
\begin{bmatrix}
a_{0,1} & a_{1,1} & \cdots & a_{N-1,1} \\
a_{0,2} & a_{1,2} & \cdots & a_{N-1,2} \\
\vdots & \vdots & \ddots & \vdots \\
a_{0,\ell} & a_{1,\ell} & \cdots & a_{N-1,\ell}
\end{bmatrix}
\f]

Or, if we define as shorthand \f$ \limbd{j}{a} := \left(a_{0, j}, \dots , a_{N-1, j}\right) \f$
(this notation will be very useful in the explanations of GLWE, GLWEGadget and GGSW), we can see it as a vector:
\f[
\begin{bmatrix}
\limbd{1}{a} & \limbd{2}{a} & \dots & \limbd{\ell}{a}
\end{bmatrix}
\f]

Where every \f$ \limbd{j}{a} \in \texttt{int64_t}_N[X] \f$ is a limb, that is, a univariate polynomial.

This representation of bivariate polynomials has what is called the "prefix property":
if we omit later limbs, the result is still a representation of the same polynomial, just with less precision.
