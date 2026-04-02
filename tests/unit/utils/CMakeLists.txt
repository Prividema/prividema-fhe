
find_package(PkgConfig REQUIRED)

# Creates imported target: PkgConfig::CRITERION
pkg_check_modules(CRITERION REQUIRED IMPORTED_TARGET criterion)

set(UTEST_UTILS_SRCS
    src/ututils.c
)

add_library(PVDA_UTEST_UTILS SHARED ${UTEST_UTILS_SRCS})

# Set the output library name explicitly
set_target_properties(PVDA_UTEST_UTILS PROPERTIES OUTPUT_NAME pvdautestutils)

# Link dependencies: bivGGSW depends on GLWE, libspqlios, and COMMON
target_link_libraries(PVDA_UTEST_UTILS PRIVATE PkgConfig::CRITERION)

# Include directories
target_include_directories(PVDA_UTEST_UTILS
    PUBLIC
    include/
)
