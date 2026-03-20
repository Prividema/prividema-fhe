# bivGLWE representation

We represent a bivGLWE as a vector arranged by precision, ie degree in Y. Let's take a set of parameters : N the degree of the cyclotomic polynomial, Kappa, ie a base-2Kappa Bg = 2^{-kappa} and a precision l. 

For a Zn[X] polynomial P(X) = Sum{0,N-1}(a_i * X^i). The base-2Kappa decomposition gives :

Dec(P(X)) = limb_1(P) * Bg + limb_2(P) * Bg^2 + ... + limb_l(P) * Bg^l

where, for i in [1, l] : 

limb_i(P) = Sum{0, N-1}(  )