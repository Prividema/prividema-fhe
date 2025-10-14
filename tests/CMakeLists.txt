cmake_minimum_required(VERSION 3.10)

# Collect all test source files
file(GLOB TEST_FILES backend/*.c core/*.c schemes/*.c)

# Create test executables for each test file
foreach(test_file ${TEST_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library
    target_link_libraries(${test_name} CORE)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
endforeach()

add_executable(conv_prod Convolution_Product/conv_prod.c )

target_link_libraries(conv_prod PUBLIC CONV)
target_link_libraries(conv_prod PUBLIC libspqlios)
target_include_directories(conv_prod PRIVATE ${PROJECT_SOURCE_DIR}/build/_deps/sp-src/spqlios/arithmetic)