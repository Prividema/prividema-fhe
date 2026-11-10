# ----------------------------------------
# ADAPTER: links between prividema and the backends
#
# Converts prividema data (pvda_types.h) into backend input data,
# casts the opaque BE_* handles to the true backend types, and
# dispatches every call to the selected backend.
# ----------------------------------------

add_subdirectory(src)

if(BUILD_TESTS)
  add_subdirectory(tests)
endif()
