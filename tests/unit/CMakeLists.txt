## Collect all test source files
file(GLOB TEST_BACKEND_FILES backend/*.c)
file(GLOB TEST_COMMON_FILES common/*.c)
file(GLOB TEST_CORE_FILES core/ggsw/*.c core/glwe/*.c)
file(GLOB TEST_SCHEMES_FILES schemes/*.c)

# Function that build the test files
function(BuildUnitTests testfiles libname)
    foreach(test_file IN LISTS ${testfiles})

       message(STATUS "hey 2 hey")
        # Get filename without extension
        get_filename_component(test_name ${test_file} NAME_WE)

        # Create executable
        add_executable(${test_name} ${test_file})

        # Link libraries
        target_link_libraries(${test_name} ${libname} criterion)

        # Add include path for headers
        target_include_directories(${test_name} PRIVATE ${PROJECT_SOURCE_DIR}/include)
        
        add_test(NAME ${test_name} COMMAND --exe $<TARGET_FILE:${test_name}>)
    endforeach()
endfunction()

# Build all the test files
BuildUnitTests("${TEST_BACKEND_FILES}" BACKEND)
BuildUnitTests("${TEST_COMMON_FILES}" COMMON)
BuildUnitTests("${TEST_CORE_FILES}" CORE)
BuildUnitTests("${TEST_SCHEMES_FILES}" SCHEMES)