# prividema-fhe

Prividema-fhe is a cryptographic library designed to unify multiple homomorphic encryption schemes under a single framework.

## Requirements

- CMake version 3.10
- [spqlios-arithmetic](https://github.com/tfhe/spqlios-arithmetic) (This will be installed automatically with CMake)
- [Criterion](https://github.com/Snaipe/Criterion) if you want to if you want to run the unit tests. In this case make sur that pkg-config is also installed for portability between different OS.
- For Windows users, also make sure that [Cryptography API: Next Generation](https://learn.microsoft.com/en-us/windows/win32/seccng/cng-portal) is installed, as it is required for the random number generator.

## To compile and run

Here's an example of how to build and run the test in `core`, from the root of the project :

```bash
# Legacy method
mkdir build;
cd build;
cmake .. -DBUILD_TESTS=ON;
make;
core/tests/hello
```

```bash
# Modern method 
cmake -S . -B build -DBUILD_TESTS=ON;
cmake --build build;
build/core/tests/hello
```

### CMake options :

- `BUILD_TESTS`: Build the test files.
- `BUILD_GGSW` : Build the ggsw library inside libcore.
- `ENABLE_DEBUG` : Enables additional debug prints.
- `BUILD_DOCS` : Build the Documentation.

## Structure 

## Security

## Docker

## References
