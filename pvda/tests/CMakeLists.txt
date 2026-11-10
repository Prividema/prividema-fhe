# ------------------------------------------------------------
# PVDA unit tests
# ------------------------------------------------------------

file(GLOB TEST_COMMON_FILES  common/*.c)
file(GLOB TEST_GGSW_FILES    core/ggsw/*.c)
file(GLOB TEST_GLWE_FILES    core/glwe/*.c)
file(GLOB TEST_SCHEMES_FILES schemes/*.c)

add_subdirectory(utils)

BuildUnitTests(TEST_COMMON_FILES  PVDA_COMMON)
BuildUnitTests(TEST_GGSW_FILES    PVDA_CORE)
BuildUnitTests(TEST_GLWE_FILES    PVDA_CORE)
BuildUnitTests(TEST_SCHEMES_FILES PVDA_SCHEMES)
