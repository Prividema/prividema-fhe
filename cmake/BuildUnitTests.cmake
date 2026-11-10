# ------------------------------------------------------------
# Unit tests helpers, shared by pvda/tests and adapter/tests
# ------------------------------------------------------------

# Criterion (portable, pkg-config only)
find_package(PkgConfig REQUIRED)

# Creates imported target: PkgConfig::CRITERION
# GLOBAL so that every test directory can link to it.
pkg_check_modules(CRITERION REQUIRED IMPORTED_TARGET GLOBAL criterion)

if(APPLE)
  set(PVDA_TESTS_NEED_MATH_LIB FALSE)
else()
  set(PVDA_TESTS_NEED_MATH_LIB TRUE)
endif()

# ------------------------------------------------------------
# BuildUnitTests(<list variable of test files> <library to test>)
#
# Builds one executable per test file, linked to the tested library,
# the test utilities and Criterion, and registers it with CTest.
# ------------------------------------------------------------
function(BuildUnitTests testfiles libname)

  foreach(test_file IN LISTS ${testfiles})
    get_filename_component(test_name ${test_file} NAME_WE)

    add_executable(${test_name} ${test_file})

    # Project includes only
    target_include_directories(${test_name} PRIVATE
            ${PVDA_INCLUDE_DIR}
            ${ADAPTER_INCLUDE_DIR}
        )

    # Link libraries
    target_link_libraries(${test_name} PRIVATE
            ${libname}
            PVDA_UTEST_UTILS
            PVDA_CORE
            # Criterion via pkg-config
            PkgConfig::CRITERION
        )

    # Math lib (Linux only)
    if(PVDA_TESTS_NEED_MATH_LIB)
      target_link_libraries(${test_name} PRIVATE m)
    endif()

    # Register test
    add_test(
            NAME ${test_name}
            COMMAND $<TARGET_FILE:${test_name}>
        )

  endforeach()

endfunction()
