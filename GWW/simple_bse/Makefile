# Makefile for simple_bse
#
include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
	$(MOD_FLAG)../../GWW/gww

#<<<<<<< HEAD
##location of needed libraries
#LIBOBJS= ../../iotk/src/libiotk.a  \
#	../../clib/clib.a
#
#=======
#>>>>>>> 40cae61e51e5d23399a72a97ff6a42645c3af4a3
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

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a \
	../../KS_Solvers/libks_solvers.a \
        ../../LAXlib/libqela.a ../../GScratch/libgscratch.a ../../UtilXlib/libutil.a

#<<<<<<< HEAD
#TLDEPS=bindir mods libs libiotk
#=======
TLDEPS=bindir mods libs gwwlib
#>>>>>>> 40cae61e51e5d23399a72a97ff6a42645c3af4a3

all : tldeps simple_bse.x  libsimple_exc.a

simple_bse.x : simple_bse.o $(EXCOBJS) $(LIBOBJS) $(QEMODS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		simple_bse.o $(EXCOBJS) $(QEMODS) $(LIBOBJS)  $(QELIBS)
	- ( cd ../../bin; ln -fs ../GWW/simple_bse/$@ . )


libsimple_exc.a : $(EXCOBJS)
	$(AR) $(ARFLAGS) $@ $(EXCOBJS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f -v simple_bse.x  *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
