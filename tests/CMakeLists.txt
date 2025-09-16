cmake_minimum_required(VERSION 3.10)

# Collect all test source files
file(GLOB TEST_FILES *.c)

# Create test executables for each test file
foreach(test_file ${TEST_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library
    target_link_libraries(${test_name} core)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
endforeach()
