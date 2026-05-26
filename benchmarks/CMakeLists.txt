cmake_minimum_required(VERSION 3.10)

set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable Google Benchmark's internal tests" FORCE)
# Google bencmakrk
include(FetchContent)
FetchContent_Declare(
  benchmark
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        v1.9.5
  CMAKE_ARGS     -DBENCHMARK_ENABLE_GTEST_TESTS=OFF -DBENCHMARK_ENABLE_ASSEMBLY_TESTS=OFF
)
FetchContent_MakeAvailable(benchmark)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)

# Collect all test source files
file(GLOB BENCH_FILES *.c *.cpp)

# Create test executables for each test file
foreach(bench_file ${BENCH_FILES})
  # Get filename without extension
  get_filename_component(bench_name ${bench_file} NAME_WE)

  # Create executable
  add_executable(${bench_name} ${bench_file})

  # Link with the main project library
  target_link_libraries(${bench_name} PRIVATE PVDA_COMMON PVDA_GLWE PVDA_GGSW PVDA_BACKEND libspqlios benchmark::benchmark_main)

  # Add include path for headers
  target_include_directories(${bench_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)

endforeach()

add_executable(bench_all ${BENCH_FILES})
target_link_libraries(bench_all PRIVATE PVDA_COMMON PVDA_GLWE PVDA_GGSW PVDA_BACKEND libspqlios benchmark::benchmark_main)
target_include_directories(bench_all PRIVATE ${PROJECT_SOURCE_DIR}/include)
