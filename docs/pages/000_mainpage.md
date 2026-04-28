# Prividema-fhe documentation {#mainpage}

Prividema-fhe is a cryptographic library designed to unify multiple homomorphic encryption schemes under a single framework.
When complete, it will contain:

- A modular backend system that can efficiently delegate raw computations to:
  - CPU libraries like SPQLIOS, that, for example, use AVX extensions
  - GPU libraries
  - FPGAs
  - ...
- A core implementing bivariate (base-2K) polynomial arithmetic and
  implementations of GLWE, GLWEGadget and GGSW using said arithmetic
- A scheme layer implementing:
  - BFV/BGV
  - TFHE
  - CKKS
  - Scheme-switching between the above (CHIMERA)
  - Encoding and decoding functions from raw data to/from the polynomial representations that the schemes use

## Library Structure

The library is to be divided in the following layers:

- Backend : will contain an abstraction layer over the underlying library or hardware that is used for heavy optimisation.
  At present only spqlios is supported
- Common: Utility code, functions that belong to no particular scheme/problem/FHE concept.
-Core: Where the code for basic mathematical constructs will go
  - GLWE: functions and code for GLWE operations
  - GGSW: functions and code for GGSW and (related) GLWEGadget opeartions
- Schemes: The different FHE schemes that can be implemented using the above problems
  
Its implementation in C allows for close-to-the metal optimisations
and maximum portability, due to the spread of C toolchains as well
as ability to call C functions (via an FFI) from most other programming languages.
Additionally, the library is structured in different layers that can be imported
independently, providing the developer a choice on the level of abstraction
that their application requires.

\image latex block.pdf "Block representation of the library's layers" width=0.55\textwidth

\image html block.svg "Block representation of the library's layers"

## Building

Instructions for building the library can be found in @ref building_lib

## Testing and benchmarks

The library currently contains both a test suite and some benchmarks.
The tests are implemented using the Criterion C framework and most of them have parametrized problem parameters (\N, \K, etc.).
The benchmarks require a relatively modern C++ compiler as they use Google's benchmark library.

The instructions can be found in @ref build_tests
