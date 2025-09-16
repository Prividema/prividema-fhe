cmake_minimum_required(VERSION 3.10)
project(prividema-fhe C)

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Add include directory
include_directories(include)

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)

