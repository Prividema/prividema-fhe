# prividema-fhe

Prividema-fhe is a cryptographic library designed to unify multiple homomorphic encryption schemes under a single framework.

## Requirements

- CMake version 3.10
- [spqlios-arithmetic](https://github.com/tfhe/spqlios-arithmetic) (This will be installed automatically with CMake)
- [Criterion](https://github.com/Snaipe/Criterion) if you want to if you want to run the unit tests. In this case make sure that pkg-config is also installed for portability between different OS.
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

> [!WARNING]
> Current Criterion versions have a bug that has tests disappear in Release builds.
> See [Issue #43](https://github.com/Prividema/prividema-fhe/issues/43) and upstream [Criterion #590](https://github.com/Snaipe/Criterion/issues/590).
> [Criterion #588](https://github.com/Snaipe/Criterion/pull/588) is a fix PR
> that was filed some months before we encountered this issue, but to this day (2026-04-20) it has not been merged.
> You should either manually apply the patch in the PR to your Criterion installation (recommended if you plan on developing the library),
> or, failing that, run tests only on builds compiled in CMake Debut release type
> (use -DCMAKE_BUILD_TYPE=Debug when invoking CMake for configuration).

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

#### Running the tests

It might be helpful to know that

```bash
ctest
```

Inside the build directory will automatically run all tests.
Make sure that you have abided by the above warning, as ctest will not by default
tell you how many tests have been run in each file (and you might thing that everything
passes when in reality no tests are being run at all).

## Structure

The library is to be divided in the following layers:

- Backend (TODO): will contain an abstraction layer over the underlying library or hardware that is used for heavy optimisations
- Common: Utility code, functions that belong to no particular scheme/problem/FHE concept.
-Core: Where the code for basic mathematical constructs will go
  - GLWE: functions and code for GLWE operations
  - GGSW: functions and code for GGSW and (related) GLWEGadget opeartions
- Schemes (TODO): The different FHE schemes that can be implemented using the above problems
  
## Security

## Docker

TODO/TBD

## References
