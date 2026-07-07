############################################################
# This is based on experience from thea Odyssey cluster
# at the University of Tokyo.
############################################################

# Stop if OFFLOAD enabled, as this is not tested
if(QE_ENABLE_OFFLOAD)
	message(FATAL_ERROR "The current setup for the Fujitsu fortran compiler does NOT "
		            "support the use of offload devices.")
endif()

############################################################
# Needed to compile EPW/src/utilities/low_lvl.f90 
# Might be useful in the future for compiler specific workarounds.
############################################################
qe_add_global_compile_definitions(__FUJITSU)

############################################################
# Set flag to allow allocating assignment,
# and -Kfast for code optiomizations.
# this is not default for the Fujitsu compiler.
############################################################

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Nalloc_assign -Kfast")

############################################################
# Suggest BLAS and LAPACK libraries
############################################################

if(QE_ENABLE_OPENMP AND NOT BLA_VENDOR)
   set(BLA_VENDOR Fujitsu_SSL2BLAMPSVE)
elseif(NOT BLA_VENDOR)
   set(BLA_VENDOR Fujitsu_SSL2SVE)
endif()


############################################################
# Note, that for SCALAPACK support to work on Odyssey
# it requires setting both QE_ENABLE_SCALAPACK
# and QE_ENABLE_OPENMP as the Fujitsu SCALAPACK
# library has only been compiled together with openMP
############################################################
