add_library(SCHEME scheme.c)
set_target_properties(SCHEME PROPERTIES OUTPUT_NAME scheme)
target_include_directories(SCHEME INTERFACE "${PROJECT_SOURCE_DIR}/include/scheme")
