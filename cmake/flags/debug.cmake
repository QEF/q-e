# Provide a human-friendly debug experience by keeping symbols and stack frames

set(flags "-Og -g -fno-omit-frame-pointer -fno-optimize-sibling-calls")

# Debug
set(CMAKE_C_FLAGS_DEBUG   "${flags}" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_DEBUG "${flags}" CACHE STRING "")

unset(flags)
