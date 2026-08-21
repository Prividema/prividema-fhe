# Prividema-fhe documentation {#mainpage}

Prividema-fhe is a cryptographic library designed to unify multiple homomorphic encryption schemes under a single framework.
When complete, it will contain:

- A modular backend system that can efficiently delegate raw computations to:
  - CPU libraries like [spqlios-arithmetic](https://github.com/tfhe/spqlios-arithmetic), that, for example, use AVX extensions
  - GPU libraries
  - FPGAs
  - ...
- A core implementing bivariate (base-2K) polynomial arithmetic and
  implementations of GLWE, GLWEGadget and GGSW using said arithmetic
- A scheme layer implementing:
  - HEInteger (HEInt)
  - HEBits  
  - HEFixedPoint (HEFP)
  - Scheme-switching between the above (CHIMERA)
  - Encoding and decoding functions from raw data to/from the polynomial representations that the schemes use

## Library Structure

The prividema-lib library is designed as a modular and layered framework to support a wide range of privacy-preserving applications. Some applications may require direct access to low-level mathematical primitives such as bivariate polynomial operations, GLWE, or GGSW, while others may benefit from higher-level implementations of complete FHE scheme abstractions  over fixed point values, integer values or circuit/LUT-based operations. This layered design allows developers to use only the components required for their specific use case, reducing complexity and improving flexibility.

The library is organized into the following main layers:

- Backend: This layer provides an abstraction over the underlying computational libraries and hardware used for performance-critical operations.
It is responsible for heavy optimizations such as FFT/NTT computations and external products.
At present, the library supports the spqlios-arithmetic spqlios-arithmetic backend, with future extensions planned for GPU- and FPGA-based acceleration from WP5.
The layer also provides random number generation, since generating it efficiently is a platform-dependent problem.
- Common: This layer contains utility code and shared functions that are not specific to any particular FHE scheme or cryptographic construction,
such as logging, error handling, and general helper functions.
- Core: This layer implements the fundamental mathematical structures and low-level cryptographic primitives:
  - GLWE: functions and data structures for GLWE ciphertexts and operations
  - GGSW / GLWEGadget: functions and implementations related to GGSW ciphertexts and GLWEGadget constructions
- Schemes: This layer contains the implementation of complete FHE schemes built on top of the core primitives, including HEFixedPoint, HEIntegers, HEBits, and scheme-switching mechanisms such as Chimera.
- Applications: This layer implements the different privacy-preserving circuits required by the project such as Private Information Retrieval (PIR) and biometric authentication.

The library is implemented in C to enable low-level optimizations (using AVX). The widespread availability of C toolchains and the possibility of calling C functions through foreign function interfaces (FFI) from most programming languages make the library suitable for integration into diverse software environments. In addition, the independent structure of the layers allows developers to choose the appropriate level of abstraction for their applications while maintaining high performance and interoperability.
The following figure shows the structure of the library.

\image latex diagram.pdf "A diagram of the library structure" width=\textwidth

\image html diagram.svg "A diagram of the library structure"

## Security

Prividema-fhe implements RLWE- and GLWE-based cryptographic primitives and follows the standard IND-CPA security model when configured with appropriate parameters.

To ensure adequate security levels, parameter selection should be performed carefully according to the target use case and security requirements. We recommend using the [Lattice Estimator](https://github.com/malb/lattice-estimator) to evaluate and validate parameter choices against known lattice attacks.

## Building

Instructions for building the library can be found in @ref building_lib

## Testing and benchmarks

The library currently contains both a test suite and some benchmarks.
The tests are implemented using the Criterion C framework and most of them have parametrized problem parameters (\N, \K, etc.).
The benchmarks require a relatively modern C++ compiler as they use Google's benchmark library.

The instructions can be found in @ref build_tests .
