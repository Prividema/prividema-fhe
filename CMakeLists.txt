cmake_minimum_required(VERSION 3.10)
project(prividema-fhe)

# ----------------------------------------
# Include FetchContent for external dependencies
# Enable testing framework
# ----------------------------------------
include(FetchContent)
include(CTest)
enable_testing()

# ----------------------------------------
# Set C standard
# ----------------------------------------
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ----------------------------------------
# Build options
# ----------------------------------------
# ----------------------------------------
# Build Options
# ----------------------------------------
# Note on CMake cache:
# When you first configure the project with specific options (e.g., -DBUILD_TESTS=ON),
# the values are stored in CMake's cache at build/CMakeCache.txt.
# If you later re-run CMake without specifying the option, it will reuse the cached value,
# which may not be what you expect.
#
# Example:
#   cmake .. -DBUILD_TESTS=ON   # BUILD_TESTS is set to ON
#   cmake ..                    # BUILD_TESTS will still be ON due to cache
#
# To reset, either:
#   - Use -U option to remove cached variables
#   - Delete the CMakeCache.txt file manually

option(BUILD_TESTS "Build the Tests" ON)
option(BUILD_EXAMPLES "Build the examples" ON)
option(BUILD_BENCHMARKS "Build the benchmarks" ON)
option(BUILD_GGSW "Build the ggsw library inside libcore" ON)
option(ENABLE_DEBUG "Enable additional prints for debugging purpose" OFF)
option(BUILD_DOCS "Build the documentation using Doxygen" OFF)

# ----------------------------------------
# Debug print configuration
# ----------------------------------------
if(ENABLE_DEBUG)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-fsanitize=address -fno-omit-frame-pointer)
    add_link_options(-fsanitize=address)
  endif()
  message(STATUS "Debug Prints Enabled")
  add_compile_definitions(DEBUG=1)
else()
  add_compile_definitions(DEBUG=0)
endif()

# ----------------------------------------
# Import spqlios library
# ----------------------------------------
FetchContent_Declare(
    spqlios
    GIT_REPOSITORY https://github.com/tfhe/spqlios-arithmetic.git
    GIT_TAG main
)
# Disable spqlios internal tests
set(ENABLE_TESTING OFF CACHE BOOL "Disable Tests" FORCE)
FetchContent_MakeAvailable(spqlios)

# ----------------------------------------
# Add main source directory
# ----------------------------------------
add_subdirectory(src)

# ----------------------------------------
# Criterion installation instructions
# ----------------------------------------

# UBUNTU USERS
# ------------------------
# 1. Install Meson (required for building Criterion):
#      sudo apt install meson
# 2. Clone and build Criterion (optional if you want latest version):
#      git clone https://github.com/Snaipe/Criterion.git
#      cd Criterion/
#      mkdir build
#      meson setup build
#      sudo meson install -C build
# 3. Verify installation:
#      ls /usr/include/criterion
#      ls /usr/lib/libcriterion*
# Notes:
#  - Linking math library (-lm) may be required on Linux

# MAC USERS
# ------------------------
# 1. Install Criterion using Homebrew (recommended):
#      brew install criterion
# 2. Ensure Meson is installed (required by Criterion):
#      brew install meson
# 3. Verify installation:
#      # Apple Silicon (arm64)
#      ls /opt/homebrew/include/criterion
#      ls /opt/homebrew/lib/libcriterion*
#      # Intel macOS
#      ls /usr/local/include/criterion
#      ls /usr/local/lib/libcriterion*
# 4. Optional: build from source
#      git clone https://github.com/Snaipe/Criterion.git
#      cd Criterion/
#      mkdir build
#      meson setup build
#      meson install -C build

# ----------------------------------------
# Build unit tests if enabled
# ----------------------------------------
if(BUILD_TESTS)
  message(STATUS "Building Test Files")
  add_subdirectory(tests)
endif()

# ----------------------------------------
# Build examples if enabled
# ----------------------------------------
if(BUILD_EXAMPLES)
  message(STATUS "Building Example Files")
  add_subdirectory(examples)
endif()

if(BUILD_BENCHMARKS)
  message(STATUS "Building Benchmark Files")
  add_subdirectory(benchmarks)
endif()

# ----------------------------------------
# Doxygen documentation
# ----------------------------------------
if(BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(doc_html ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMENT "Generating HTML documentation from Markdown"
            VERBATIM
        )
  else()
    message(WARNING "Doxygen needs to be installed to generate the documentation")
  endif()
endif()
