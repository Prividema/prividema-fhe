# ------------------------------------------------------------
# Platform-specific configuration
# ------------------------------------------------------------

if(APPLE)
    # Homebrew prefix detection
    if(CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
        set(HOMEBREW_PREFIX /opt/homebrew)
    else()
        set(HOMEBREW_PREFIX /usr/local)
    endif()

    set(CRITERION_INCLUDE_DIR ${HOMEBREW_PREFIX}/include)
    set(CRITERION_LIB_DIR     ${HOMEBREW_PREFIX}/lib)

    # macOS does NOT need libm
    set(NEED_MATH_LIB FALSE)
else()
    # Linux defaults
    set(CRITERION_INCLUDE_DIR /usr/include)
    set(CRITERION_LIB_DIR     /usr/lib)

    # Linux needs libm
    set(NEED_MATH_LIB TRUE)
endif()

# ------------------------------------------------------------
# Collect all test source files
# ------------------------------------------------------------

file(GLOB TEST_BACKEND_FILES backend/*.c)
file(GLOB TEST_COMMON_FILES  common/*.c)
file(GLOB TEST_CORE_FILES
    core/ggsw/*.c
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
        #  - Criterion headers
        #  - SPQLIOS headers
        target_include_directories(${test_name} PRIVATE
            ${PROJECT_SOURCE_DIR}/include
            ${CRITERION_INCLUDE_DIR}
            ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/
        )

        # Library search path for Criterion
        target_link_directories(${test_name} PRIVATE
            ${CRITERION_LIB_DIR}
        )

        # Link project library + Criterion + SPQLIOS
        target_link_libraries(${test_name} PRIVATE
            ${libname}
            criterion
            libspqlios   # <-- Add SPQLIOS here
        )

        # Link math library only where required (Linux)
        if(NEED_MATH_LIB)
            target_link_libraries(${test_name} PRIVATE m )
        endif()

        # Ensure Criterion is found at runtime on macOS
        if(APPLE)
            set_target_properties(${test_name} PROPERTIES
                BUILD_RPATH ${CRITERION_LIB_DIR}
            )
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
BuildUnitTests(TEST_CORE_FILES    CORE)
BuildUnitTests(TEST_SCHEMES_FILES SCHEMES)
