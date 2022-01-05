# Makefile for head (head.x)
# 
# Author G. Stenuit
#
include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
	  $(MOD_FLAG)../../PW/src \
	  $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../dft-d3 \

HEADOBJS = \
bcast_gwl_input.o \
close_gwl.o \
openfilq_gwl.o \
gwl_readin.o \
lanczos_k.o \
solve_head.o 

QEMODS = $(BASEMODS)
LIBPWPH =  ../../PHonon/PH/libph.a ../../LR_Modules/liblrmod.a \
	   ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a


LIBMIN= ../minpack/minpacklib.a

TLDEPS= phlibs 

all : tldeps head.x  

head.x : head.o $(HEADOBJS) $(LIBPWPH) $(LIBMIN)
	$(LD) $(LDFLAGS) -o head.x head.o \
		 $(HEADOBJS) $(LIBPWPH)  $(QEMODS) $(LIBMIN) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../GWW/head/head.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

include make.depend
