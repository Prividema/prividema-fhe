# prividema-fhe
Prividema-fhe is a cryptographic library designed to unify multiple homomorphic encryption schemes under a single framework.

# Requirements
-spqlios-arithmetic
-cmake version 3.10

# To compile and run
Here's an example of how to build and run the test in `core`, from the root of the project :
```bash
# Old method
mkdir build
cd build
cmake .. -DBUILD_TESTS=ON
make
core/tests/hello
```

```bash
# New method 
cmake -S . -B build -DBUILD_TESTS=ON -DENABLE_TESTING=OFF
cmake --build build
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
