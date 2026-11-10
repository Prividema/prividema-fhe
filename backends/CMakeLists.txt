# ----------------------------------------
# Computation backends
#
# Every backend prividema can run on is imported here. The adapter
# (see adapter/) is the only part of the project that links to them.
# ----------------------------------------

# spqlios-arithmetic (git submodule)
# Disable spqlios internal tests
set(ENABLE_TESTING OFF CACHE BOOL "Disable Tests" FORCE)
add_subdirectory(spqlios-arithmetic)
