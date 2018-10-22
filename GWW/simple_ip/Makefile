# Makefile for simple_ip

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../pw4gww $(MOD_FLAG)../../PW/src

#<<<<<<< HEAD
##location of needed libraries
#LIBOBJS= ../../iotk/src/libiotk.a  \
#         ../../clib/clib.a

#=======
#>>>>>>> 40cae61e51e5d23399a72a97ff6a42645c3af4a3
SIMPLEIPOBJS = \
simple_ip_objects.o \
start_end.o \
input_simple_ip.o \
tetra_mod1.o \
interpolation.o \
diagonalization.o \
dielectric.o


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a \
        ../../KS_Solvers/libks_solvers.a \
	../../LAXlib/libqela.a ../../GScratch/libgscratch.a ../../UtilXlib/libutil.a
PWOBJS = ../../PW/src/libpw.a

LIBMIN=

TLDEPS=pwlibs

all : tldeps simple_ip.x  libsimple_ip.a

simple_ip.x : simple_ip.o $(SIMPLEIPOBJS) $(LIBOBJS) $(PWOBJS) $(QEMODS)  $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
		simple_ip.o $(SIMPLEIPOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBMIN) $(QELIBS)
	- ( cd ../../bin; ln -fs ../GWW/simple_ip/$@ . )

libsimple_ip.a : $(SIMPLEIPOBJS)
	$(AR) $(ARFLAGS) $@ $(SIMPLEIPOBJS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f -v simple_ip.x  *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
