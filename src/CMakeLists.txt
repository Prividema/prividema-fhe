# ----------------------------------------
# Build project libraries
# ----------------------------------------

# Backend library
add_subdirectory(backend)

# Common library
add_subdirectory(common)

# Core library (interface library containing bivGLWE and bivGGSW modules)
add_subdirectory(core)

# Schemes library
add_subdirectory(schemes)
