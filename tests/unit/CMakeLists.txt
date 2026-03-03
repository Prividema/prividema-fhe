# ------------------------------------------------------------
# Criterion (portable, pkg-config only)
# ------------------------------------------------------------

find_package(PkgConfig REQUIRED)

# Creates imported target: PkgConfig::CRITERION
pkg_check_modules(CRITERION REQUIRED IMPORTED_TARGET criterion)

# ------------------------------------------------------------
# Platform-specific configuration
# ------------------------------------------------------------

<<<<<<< HEAD
=======
# Fetch criterion with pkg-config
find_package(PkgConfig REQUIRED)

pkg_check_modules(
    CRITERION
    REQUIRED
    IMPORTED_TARGET
    criterion
)

>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
if(APPLE)
    set(NEED_MATH_LIB FALSE)
else()
    set(NEED_MATH_LIB TRUE)
endif()

# ------------------------------------------------------------
# Collect all test source files
# ------------------------------------------------------------

file(GLOB TEST_BACKEND_FILES backend/*.c)
file(GLOB TEST_COMMON_FILES  common/*.c)
<<<<<<< HEAD
file(GLOB TEST_GGSW_FILES    core/ggsw/*.c)
file(GLOB TEST_GLWE_FILES    core/glwe/*.c)
=======
file(GLOB TEST_GGSW_FILES core/ggsw/*.c)
file(GLOB TEST_GLWE_FILES core/glwe/*.c)
>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
file(GLOB TEST_SCHEMES_FILES schemes/*.c)

# ------------------------------------------------------------
# Function to build unit tests
# ------------------------------------------------------------

function(BuildUnitTests testfiles libname)

    foreach(test_file IN LISTS ${testfiles})
<<<<<<< HEAD
=======

>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
        get_filename_component(test_name ${test_file} NAME_WE)

        add_executable(${test_name} ${test_file})

<<<<<<< HEAD
=======
        # Project includes only
>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
        target_include_directories(${test_name} PRIVATE
            ${PROJECT_SOURCE_DIR}/include
            ${PROJECT_SOURCE_DIR}/build/_deps/spqlios-src/spqlios/arithmetic/
        )

<<<<<<< HEAD
=======
        # Link libraries
>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
        target_link_libraries(${test_name} PRIVATE
            ${libname}
            GLWE
            COMMON
<<<<<<< HEAD
            PkgConfig::CRITERION
            libspqlios
        )

=======
            libspqlios

            # Criterion via pkg-config
            PkgConfig::CRITERION
        )

        # Math lib (Linux only)
>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
        if(NEED_MATH_LIB)
            target_link_libraries(${test_name} PRIVATE m)
        endif()

<<<<<<< HEAD
=======
        # Register test
>>>>>>> 316cc65f17a2f8568698d49a7cf35d246ab4600d
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
BuildUnitTests(TEST_SCHEMES_FILES SCHEMES)