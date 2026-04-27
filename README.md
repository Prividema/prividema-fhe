# prividema-fhe

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

## Bivariate (base-2K) polynomial representation

Most FHE library implementations make use of Residue Number System (RNS) representation to represent
and efficiently compute the large polynomials that make up current lattice-based FHE schemes.

TODO

## Library Structure

The library is to be divided in the following layers:

- Backend (TODO): will contain an abstraction layer over the underlying library or hardware that is used for heavy optimisations
- Common: Utility code, functions that belong to no particular scheme/problem/FHE concept.
-Core: Where the code for basic mathematical constructs will go
  - GLWE: functions and code for GLWE operations
  - GGSW: functions and code for GGSW and (related) GLWEGadget opeartions
- Schemes (TODO): The different FHE schemes that can be implemented using the above problems
  
Its implementation in C allows for close-to-the metal optimisations
and maximum portability, due to the spread of C toolchains as well
as ability to call C functions (via an FFI) from most other programming languages.
Additionally, the library is structured in different layers that can be imported
independently, providing the developer a choice on the level of abstraction
that their application requires.

## Security

TODO

## Building and testing the library

Instructions for building and testing the library can be found in [BUILDING.md](BUILDING.md).

## Docker

A Docker image for building and testing the library will be provided in the future.

## References

TODO
