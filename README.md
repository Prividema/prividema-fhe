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

### CMake options

- `BUILD_TESTS`: Build the test files.
- `ENABLE_DEBUG` : Enables additional debug prints and enables sanitizers.
- `BUILD_DOCS` : Build the Documentation.
- `BUILD_EXAMPLES`: Build the example executables (TODO)
- `BUILD_BENCHMARKS`: Build the benchmarks (requires C++ compiler)
- `BUILD_NATIVE`: Enables -march=native optimisations, makes output NON-portable
- `PROFILING_OPTIONS`: Adds debugging options and info even in Release builds, intended for profiling traces

### Building the tests

#### Criterion installation instructions

##### UBUNTU USERS

 1. Install Meson (required for building Criterion):

        sudo apt install meson

 1. Clone and build Criterion (optional if you want latest version):

        git clone <https://github.com/Snaipe/Criterion.git>
        cd Criterion/
        mkdir build
        meson setup build
        sudo meson install -C build

 1. Verify installation:

        ls /usr/include/criterion
        ls /usr/lib/libcriterion*

 Notes:

- Linking math library (-lm) may be required on Linux

##### MAC USERS

 1. Install Criterion using Homebrew (recommended):

        brew install criterion

 2. Ensure Meson is installed (required by Criterion):

        brew install meson

 3. Verify installation:

     _Apple Silicon (arm64)_

        ls /opt/homebrew/include/criterion
        ls /opt/homebrew/lib/libcriterion*

    _Intel macOS_

        ls /usr/local/include/criterion
        ls /usr/local/lib/libcriterion*

 4. Optional: build from source

        git clone <https://github.com/Snaipe/Criterion.git>
        cd Criterion/
        mkdir build
        meson setup build
        meson install -C build

## Structure

## Security

## Docker

## References
