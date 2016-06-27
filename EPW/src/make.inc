#=======================================
# This file was copied and modified from the W90 directory inside QE-4.0.3
# For use inside a Quantum Espresso Dist
#=======================================

# read in the QE config
include ../../make.inc

# here we use sed to re-write the library paths. This is only relevant
# if we are using the Q-E internal lapack and blas. Q-E refers to them
# as ../flib/*.a, when what we need is ../../flib*.a

#LIBS=$(SCALAPACK_LIBS) $(shell echo '$(LAPACK_LIBS)' | sed -e 's/..\/flib/..\/..\/flib/g') $(shell echo  '$(BLAS_LIBS)' | sed -e 's/..\/flib/..\/..\/flib/g')  $(FFT_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(PGPLOT_LIBS) $(LD_LIBS)

#LIBS           = $(SCALAPACK_LIBS) $(LAPACK_LIBS) $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(MASS_LIBS) $(LD_LIBS)

# If this causes a problem, and you aren't using the Q-E lapack and blas
# uncomment this simpler version
#LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(MASS_LIBS) $(PGPLOT_LIBS) $(LD_LIBS)

#
# In the case of gfortran we need to explicitly allow for lines longer than 132 chars
# by adding the option "-ffree-form -ffree-line-length-none" in F90FLAGS
# FG - June 2014
#
ifneq ($(filter -D__GFORTRAN ,$(DFLAGS)),)
    MPIF90 := $(MPIF90) -ffree-form -ffree-line-length-none
endif
#
# In the case of ifort we need to make sure that temp arrays don't go in the stack,
# otherwise we may have segfaults with no apparent reason (the actual reason being
# that some versions of ifort try to allocate temp arrays in a stack which is already full)
# FG - July 2014
#
COMPILER = $(shell $(MPIF90) --version | grep ifort | sed 's/ .*//' )
ifeq "$(COMPILER)" "ifort"
    MPIF90 := $(MPIF90) -heap-arrays
endif

