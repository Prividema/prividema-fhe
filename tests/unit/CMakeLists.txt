## Collect all test source files
file(GLOB TEST_BACKEND_FILES backend/*.c)
file(GLOB TEST_COMMON_FILES common/*.c)
file(GLOB TEST_CORE_FILES core/ggsw/*.c core/glwe/*.c)
file(GLOB TEST_SCHEMES_FILES schemes/*.c)

# Create Unit test executables for each test file
foreach(test_file ${TEST_BACKEND_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library and Criterion 
    target_link_libraries(${test_name} BACKEND criterion gsl)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
    add_test(NAME ${test_name} COMMAND --exe $<TARGET_FILE:${test_name}>)
endforeach()

# Create Unit test executables for each test file
foreach(test_file ${TEST_COMMON_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library and Criterion 
    target_link_libraries(${test_name} COMMON criterion gsl)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
    add_test(NAME ${test_name} COMMAND --exe $<TARGET_FILE:${test_name}>)
endforeach()

# Create Unit test executables for each test file
foreach(test_file ${TEST_CORE_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library and Criterion 
    target_link_libraries(${test_name} CORE criterion gsl)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
    add_test(NAME ${test_name} COMMAND --exe $<TARGET_FILE:${test_name}>)
endforeach()

# Create Unit test executables for each test file
foreach(test_file ${TEST_SCHEMES_FILES})
    # Get filename without extension
    get_filename_component(test_name ${test_file} NAME_WE)

    # Create executable
    add_executable(${test_name} ${test_file})

    # Link with the main project library and Criterion 
    target_link_libraries(${test_name} SCHEMES criterion gsl)

    # Add include path for headers
    target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
    
    add_test(NAME ${test_name} COMMAND --exe $<TARGET_FILE:${test_name}>)
endforeach()