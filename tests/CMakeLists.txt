cmake_minimum_required(VERSION 3.10)

## Collect all test source files
file(GLOB TEST_BACKEND_FILES backend/*.c)

# Create test executables for each test file
foreach(test_file ${TEST_BACKEND_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Include Criterion headers
    target_include_directories(${test_name} PRIVATE /usr/include)

    # Link with the main project library and Criterion
    target_link_libraries(${test_name} BACKEND criterion)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
endforeach()


## Collect all test source files
file(GLOB TEST_CORE_FILES core/*.c)

# Create test executables for each test file
foreach(test_file ${TEST_CORE_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Add directory of Criterion 
    target_include_directories(${test_name} PRIVATE /usr/include)

    # Link with the main project library and Criterion
    target_link_libraries(${test_name} CORE criterion)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
    add_test(NAME HelloTests COMMAND ${test_name})
endforeach()


## Collect all test source files
file(GLOB TEST_SCHEME_FILES schemes/*.c)

# Create test executables for each test file
foreach(test_file ${TEST_SCHEME_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library
    target_link_libraries(${test_name} SCHEME)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
endforeach()