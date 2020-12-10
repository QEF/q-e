# Auto-detect uarch of the build machine and automatically codegen for it
# Be aware that resulting binaries will be unportable even across same arch but different uarch

set(flags "-O3 -DNDEBUG -march=native -mtune=native")

# Release
set(CMAKE_C_FLAGS_RELEASE   "${flags}" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "${flags}" CACHE STRING "")

# RelWithDebInfo
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "${flags} -g" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${flags} -g" CACHE STRING "")

unset(flags)
