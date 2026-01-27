# ------------------------------------------------------------
# Platform-specific configuration
# ------------------------------------------------------------

# Fetch criterion
find_package(PkgConfig REQUIRED)
pkg_check_modules(CRITERION REQUIRED criterion)

if(APPLE)
    # macOS does NOT need libm
    set(NEED_MATH_LIB FALSE)
else()
    set(NEED_MATH_LIB TRUE)
endif()

# ------------------------------------------------------------
# Collect all test source files
# ------------------------------------------------------------

file(GLOB TEST_BACKEND_FILES backend/*.c)
file(GLOB TEST_COMMON_FILES  common/*.c)
file(GLOB TEST_GGSW_FILES
    core/ggsw/*.c
)
file(GLOB TEST_GLWE_FILES
    core/glwe/*.c
)
file(GLOB TEST_SCHEMES_FILES schemes/*.c)

# ------------------------------------------------------------
# Function to build unit tests
# ------------------------------------------------------------

function(BuildUnitTests testfiles libname)
    foreach(test_file IN LISTS ${testfiles})
        # Get test name without extension
        get_filename_component(test_name ${test_file} NAME_WE)

        # Create test executable
        add_executable(${test_name} ${test_file})

        # Include paths:
        #  - Project headers
        #  - SPQLIOS headers
        target_include_directories(${test_name} PRIVATE
            ${PROJECT_SOURCE_DIR}/include
            ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/
        )
        # Link project library + Criterion + SPQLIOS
        target_link_libraries(${test_name} PRIVATE
            ${libname}
            GLWE
            COMMON
            criterion
            libspqlios   # <-- Add SPQLIOS here
        )

        # Link math library only where required (Linux)
        if(NEED_MATH_LIB)
            target_link_libraries(${test_name} PRIVATE m )
        endif()

        # Register test with CTest
        add_test(
            NAME ${test_name}
            COMMAND $<TARGET_FILE:${test_name}>
        )
    endforeach()
endfunction()

# ------------------------------------------------------------
# Build all unit tests
# ------------------------------------------------------------

BuildUnitTests(TEST_BACKEND_FILES BACKEND)
BuildUnitTests(TEST_COMMON_FILES  COMMON)
BuildUnitTests(TEST_GGSW_FILES    CORE)
BuildUnitTests(TEST_GLWE_FILES    CORE)
BuildUnitTests(TEST_SCHEMES_FILES  SCHEMES)