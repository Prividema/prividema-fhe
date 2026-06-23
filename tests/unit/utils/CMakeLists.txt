
find_package(PkgConfig REQUIRED)

# Creates imported target: PkgConfig::CRITERION
pkg_check_modules(CRITERION REQUIRED IMPORTED_TARGET criterion)

set(UTEST_UTILS_SRCS
    src/test_utils.c
    src/stat_utils.c
    src/ggsw_utils.c
    src/glwegadget_utils.c
)

add_library(PVDA_UTEST_UTILS ${UTEST_UTILS_SRCS})

# Set the output library name explicitly
set_target_properties(PVDA_UTEST_UTILS PROPERTIES OUTPUT_NAME pvdautestutils)

target_link_libraries(PVDA_UTEST_UTILS PRIVATE PkgConfig::CRITERION PVDA_GLWE PVDA_COMMON PVDA_GGSW PVDA_BACKEND)

# Include directories
target_include_directories(PVDA_UTEST_UTILS
    PUBLIC
    include/
)
