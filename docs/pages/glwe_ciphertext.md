# Bivariate GLWE

This file documents the internal representation of a bivariate GLWE ciphertext used in the library.

## Notation

See [the bivariate polynomial documentation's memory representation section](\ref bivariate-notation) for the notation used throughout.

## Memory layout

A torus GLWE encryption of m is \f$ \glwe(m) = \left( a_0, a_1, \dots, a_{k-1}, b \right) \f$, where  \f$ a_i \in \tnx \f$.
It is common to define \f$ a_k := b \f$ for convenience.

In the library, GLWEs are laid out in memory as to also have a prefix property.
Thus, if the decomposition in limbs of eg. \f$ a_0 \f$ in limbs \f$ \limbd{j}{a_0} \f$
 is \f$ {a_0} := \left(\limbd{1}{a_0}, \dots , \limbd{\ell_A}{a_0}\right) \f$,
then \f$ \glwe(m) \f$ is encoded when \f$ l_A = l_b \f$ as the follow matrix flattened in a row-major order:

\f[
\begin{split}
\glwe(m) & = \begin{bmatrix}
  \limbd{1}{a_0} & \limbd{1}{a_1} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} \\
  \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \limbd{2}{a_{k-1}} & \limbd{2}{b} \\
  \vdots        & \vdots        & \ddots & \vdots & \vdots \\
  \limbd{l_A-1}{a_0} & \limbd{l_A-1}{a_1} & \dots & \limbd{l_A-1}{a_{k-1}} & \limbd{l_b-1}{b} \\
  \limbd{l_A}{a_0} & \limbd{l_A}{a_1} & \dots & \limbd{l_A}{a_{k-1}} & \limbd{l_b}{b}
\end{bmatrix} \\
&  \\
 & = \begin{bmatrix}
  \limbd{1}{a_0} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} & \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \dots & \limbd{l_A}{a_{k-1}} & \limbd{l_b}{b}
\end{bmatrix}
\end{split}
\f]

Prividema-lib also supports \f$ l_A = l_b + 1\f$ (only other possibility), in which case, the memory layout is:

\f[

\begin{split}
\glwe(m) & = \begin{bmatrix}
  \limbd{1}{a_0} & \limbd{1}{a_1} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} \\
  \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \limbd{2}{a_{k-1}} & \limbd{2}{b} \\
  \vdots        & \vdots        & \ddots & \vdots & \vdots \\
  \limbd{l_A-1}{a_0} & \limbd{l_A-1}{a_1} & \dots & \limbd{l_A-1}{a_{k-1}} & \limbd{l_b}{b} \\
  \limbd{l_A}{a_0} & \limbd{l_A}{a_1} & \dots & \limbd{l_A}{a_{k-1}} &
\end{bmatrix}\\
& \\
& =
  \begin{bmatrix}
  \limbd{1}{a_0} & \limbd{1}{a_1} & \dots & \limbd{1}{a_{k-1}} & \limbd{1}{b} &
  \limbd{2}{a_0} & \limbd{2}{a_1} & \dots & \dots & \limbd{l_A}{a_{k-1}}
\end{bmatrix}
\end{split}
\f]
