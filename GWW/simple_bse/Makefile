# Makefile for simple_bse
#
include ../../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
	$(MOD_FLAG)../../GWW/gww 

EXCOBJS = \
simple_object.o \
start_end.o \
input_simple_exc.o \
simple_eigen.o \
diago_exc_sd.o \
hamiltonian.o \
derived_object.o \
diago_exc_cg.o \
spectrum.o \
lanczos.o \
build_eemat.o

QEMODS = $(BASEMODS)

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS= mods gwwlib
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../../install/tldeps

all : simple_bse.x  libsimple_exc.a

simple_bse.x : simple_bse.o $(EXCOBJS) $(QEMODS)  
	$(MPIF90) $(LDFLAGS) -o $@ \
		simple_bse.o $(EXCOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../GWW/simple_bse/$@ . )


libsimple_exc.a : $(EXCOBJS)
	$(AR) $(ARFLAGS) $@ $(EXCOBJS)

clean :
	- /bin/rm -f -v simple_bse.x  *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
