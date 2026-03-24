cmake_minimum_required(VERSION 3.10)

# Collect all test source files
file(GLOB BENCH_FILES *.c)

# Create test executables for each test file
foreach(bench_file ${BENCH_FILES})
  # Get filename without extension
  get_filename_component(bench_name ${bench_file} NAME_WE)

  # Create executable
  add_executable(${bench_name} ${bench_file})

  # Link with the main project library
  target_link_libraries(${bench_name} CORE)

  # Add include path for headers
  target_include_directories(${bench_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)

endforeach()
