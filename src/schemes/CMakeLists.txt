add_library(SCHEME scheme.c)
target_include_directories(SCHEME INTERFACE "${PROJECT_SOURCE_DIR}/include/scheme")