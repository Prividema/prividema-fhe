# Prividema-fhe library structure

The Prividema-fhe library is designed to be modular to allow it to be
used for a variety of applications.
Some applications might use the mathematical primitives (GLWE, GGSW, or even
the lower level bivariate polynomial operations directly), while others
might benefit from not having to implement common schemes like TFHE, BFV, etc.

The library is therefore structured in a layered fashion, and application
developers can choose to link only with the parts that they need in their application.
The following figure shows the structure of the library.

\image latex diagram.pdf "A diagram of the library structure" width=\textwidth

\image html diagram.svg "A diagram of the library structure"
