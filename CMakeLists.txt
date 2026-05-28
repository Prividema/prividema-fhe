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
option(BUILD_BENCHMARKS "Build the benchmarks" OFF)
option(PROFILING_OPTIONS "Add options to allow for proper profiling even in realease-like builds" OFF)
option(BUILD_NATIVE "Enable march=native optimisations" OFF)
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

if(PROFILING_OPTIONS)
  add_compile_options(-fno-omit-frame-pointer -g)
endif()

if(BUILD_NATIVE)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|IntelLLVM")
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag("-march=native" COMPILER_SUPPORTS_MARCH_NATIVE)

    if(COMPILER_SUPPORTS_MARCH_NATIVE)
      add_compile_options(-march=native)
      message(STATUS "Native optimization (-march=native) enabled.")
    else()
      message(WARNING "-march=native is not supported by this compiler.")
    endif()

  elseif(MSVC)
    message(STATUS "MSVC detected. -march=native skipped (manual /arch flag required).")
  endif()

endif()

include(CheckIPOSupported)
check_ipo_supported(RESULT lto_supported OUTPUT lto_error)

if(lto_supported)
  message(STATUS "LTO/IPO enabled for optimized builds")
  # This sets the global default for all targets defined after this point
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
  add_compile_options(-flto)
  add_link_options(-flto)
else()
  message(STATUS "LTO not supported: ${lto_error}")
endif()
# ----------------------------------------
# Import spqlios library
#
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
# Build unit tests if enabled
# ----------------------------------------
if(BUILD_TESTS)
  message(STATUS "Building Test Files")
  if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message(WARNING
    "Upstream Criterion 2.4.3 (and possibly older versions) have a bug where Release (optimized) builds"
    "can have their parameterized tests be optimized away (deleted).\n"
    "PATCH YOUR LOCAL CRITERION INSTALLATION OR COMPILE IN Debug MODE IF YOU WANT TO RUN THE TESTS!"
    )
  endif()
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
