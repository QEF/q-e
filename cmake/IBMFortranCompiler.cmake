# set preprocessor specific flag
set(Fortran_PREPROCESSOR_FLAGS "-qsuffix=cpp=f90")

# set the flag to ensure symbols with underscore added for linking Fortran libraries built by other compilers.
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qextname")
