# Makefile for Modules

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src \
	  $(MOD_FLAG)../Modules \
	  $(MOD_FLAG)../EE $(MOD_FLAG)../PW $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../iotk/src/libiotk.a ../flib/flib.a \
         ../clib/clib.a ../flib/ptools.a ../multigrid/mglib.a

NEBOBJS = \
string_methods.o \
compute_scf.o \
compute_fes_grads.o 

NEBLIBS = \
input.o \
path_base.o \
path_formats.o \
path_io_routines.o \
path_opt_routines.o \
path_reparametrisation.o \
path_variables.o

QEMODS=../Modules/libqemod.a
PWOBJS= ../PW/libpw.a
EEOBJS=../EE/libee.a

TLDEPS=bindir mods libs liblapack libblas pw eelib

all : tldeps neb.x

neb.x : $(NEBOBJS) libneb.a $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	$(NEBOBJS) libneb.a $(EEOBJS) $(LIBOBJS) $(QEMODS) $(PWOBJS) $(LIBS)
	- ( cd ../bin; ln -fs ../NEB/$@ . )

libneb.a : $(NEBLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f pw.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

include make.depend
