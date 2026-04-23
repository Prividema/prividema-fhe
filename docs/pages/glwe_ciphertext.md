# Bivariate GLWE, GLWEGadget and GGSW

This file documents the internal representation of a bivariate GLWE ciphertext used in the library.

## Notation

See [the bivariate polynomial documentation's memory representation section](\ref bivariate-notation) for the notation used throughout.

## GLWE Memory layout

A torus GLWE encryption of m is \f$ \glwe(m) = \left( a_0, a_1, \dots, a_{k-1}, b \right) \f$, where  \f$ a_i \in \tnx \f$.
It is common to define \f$ a_k := b \f$ for convenience.

In the library, GLWEs are laid out in memory as to also have a prefix property.
Thus, if the decomposition in limbs of eg. \f$ a_0 \f$ in limbs \f$ \limbd{j}{a_0} \f$
 is \f$ {a_0} := \left(\limbd{1}{a_0}, \dots , \limbd{\ell_A}{a_0}\right) \f$,
then \f$ \glwe(m) \f$ is encoded when \f$ \ell_A = \ell_b \f$ as the following matrix flattened in a row-major order:

\f[
\begin{split}
\glwe(m) & = \begin{bmatrix}
  \limbd{1}{a_0} & \limbd{1}{a_1} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} \\
  \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \limbd{2}{a_{k-1}} & \limbd{2}{b} \\
  \vdots        & \vdots        & \ddots & \vdots & \vdots \\
  \limbd{\ell_A-1}{a_0} & \limbd{\ell_A-1}{a_1} & \dots & \limbd{\ell_A-1}{a_{k-1}} & \limbd{\ell_b-1}{b} \\
  \limbd{\ell_A}{a_0} & \limbd{\ell_A}{a_1} & \dots & \limbd{\ell_A}{a_{k-1}} & \limbd{\ell_b}{b}
\end{bmatrix} \\
&  \\
 & = \flatGLWE{}
\end{split}
\f]

Prividema-lib also supports \f$ \ell_A = \ell_b + 1\f$ (only other possibility), in which case, the memory layout is:

\f[

\begin{split}
\glwe(m) & = \begin{bmatrix}
  \limbd{1}{a_0} & \limbd{1}{a_1} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} \\
  \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \limbd{2}{a_{k-1}} & \limbd{2}{b} \\
  \vdots        & \vdots        & \ddots & \vdots & \vdots \\
  \limbd{\ell_A-1}{a_0} & \limbd{\ell_A-1}{a_1} & \dots & \limbd{\ell_A-1}{a_{k-1}} & \limbd{\ell_b}{b} \\
  \limbd{\ell_A}{a_0} & \limbd{\ell_A}{a_1} & \dots & \limbd{\ell_A}{a_{k-1}} &
\end{bmatrix}\\
& \\
& =
  \begin{bmatrix}
  \limbd{1}{a_0} & \limbd{1}{a_1} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} &
  \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \dots & \limbd{\ell_A}{a_{k-1}}
\end{bmatrix}
\end{split}
\f]

## GLWEGadget

A GLWEGadget of m is a family of encryptions of m by powers of \f$ 2^{-K}\f$:
\f[
\glwegad(m) = \left[\glwe(m \cdot 2^{-K}), \glwe(m \cdot 2^{-2K}), \dots , \glwe(m \cdot 2^{-\tilde{l}K})\right]
\f]

Its memory layout is the concatenation of the GLWEs that compose it, ie (again, as a row-major matrix):

\f[
\glwegad(m) = \begin{bmatrix}
\glwe(m \cdot 2^{-K}) = \flatGLWE{1}\\
\glwe(m \cdot 2^{-2K}) = \flatGLWE{2}\\\
\vdots \\
\glwe(m \cdot 2^{-\tilde{l}K}) = \flatGLWE{\tilde{l}}
\end{bmatrix}
\f]

## GGSW

Finally a GGSW ciphertext, given its secret key \f$ s = \left(s_0, ..., s_{k-1} \right) , s \in \znx \f$,
is a set of GLWEGadgets of \f$\left(-s_0 \cdot m, \dots, -s_{k-1} \cdot m , m\right)\f$, which we lay out strided
over the GLWEGadget limbs:

\f[

\ggsw(m) = \begin{bmatrix}
\glwe(-m \cdot s_0 \cdot 2^{-K}) =& \flatGLWEW{1}{0}\\
\glwe(-m \cdot s_1 \cdot  2^{-K}) =& \flatGLWEW{1}{1}\\\
&\vdots \\
\glwe(m \cdot 2^{-K}) =& \flatGLWEW{0}{k}\\\
\glwe(-m \cdot s_0 \cdot 2^{-2K}) = &\flatGLWEW{2}{0}\\
\glwe(-m \cdot s_1 \cdot 2^{-2K}) = &\flatGLWEW{2}{1}\\\
&\vdots \\
\glwe(-m \cdot s_{k-1} \cdot 2^{-\tilde{\ell_A}K}) = &\flatGLWEW{\tilde{\ell_A}}{k-1}\\\
\glwe(m \cdot 2^{-\tilde{\ell_b}K}) = & \flatGLWEW{\tilde{\ell_b}}{k}\\\
\end{bmatrix}

\f]

In this case, we can also have \f$ \ell_A = \ell_b + 1\f$ instead of the above pictured \f$ \ell_A = \ell_b\f$, in which case the last GLWE is missing.
