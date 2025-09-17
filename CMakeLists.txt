cmake_minimum_required(VERSION 3.10)
project(prividema-fhe C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Add include directory
include_directories(include)

# Add subdirectories
add_subdirectory(src)

# Add tests subdirectory if the build option BUILD_TESTS is ON
# The option is OFF by default and can be turned ON by executing cmake with -DBUILD_TESTS=ON
option(BUILD_TESTS "Build the Tests" OFF)
if(BUILD_TESTS)
    add_subdirectory(tests)
endif(BUILD_TESTS)