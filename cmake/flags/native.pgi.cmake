# Auto-detect uarch of the build machine and automatically codegen for it
# Be aware that resulting binaries will be unportable even across same arch but different uarch

set(flags "-O3 -DNDEBUG -march=native -mtune=native")
set(fflags "-fast -Mcache_align -Mlarge_arrays -mp")

# Release
set(CMAKE_C_FLAGS_RELEASE "${flags}" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "${flags} ${fflags}" CACHE STRING "")

# RelWithDebInfo
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${flags} -g" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO "${flags} ${fflags} -g" CACHE STRING "")

unset(flags)
unset(fflags)
