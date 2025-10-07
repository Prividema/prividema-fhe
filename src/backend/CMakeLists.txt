# SPQLios
add_library(BACKEND SHARED backend.c)
target_include_directories(BACKEND INTERFACE "${PROJECT_SOURCE_DIR}/include/backend")