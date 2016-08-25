# Makefile for head (head.x)
# 
# Author G. Stenuit
#
include ../../make.inc

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG)../../LR_Modules $(MOD_FLAG)../../PW/src \
	  $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../FFTXlib $(MOD_FLAG).

HEADOBJS = \
bcast_ph_input.o \
close_phq.o \
lanczos_k.o \
openfilq.o \
phq_readin.o \
solve_head.o 

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
LIBPWPH =  ../../PHonon/PH/libph.a ../..//LR_Modules/liblrmod.a ../../PW/src/libpw.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS= bindir mods libs pw ph 

all : tldeps head.x  

head.x : head.o $(LIBOBJS) $(HEADOBJS) $(LIBPWPH) $(LIBMIN)
	$(LD) $(LDFLAGS) -o head.x head.o \
		 $(HEADOBJS) $(LIBPWPH)  $(QEMODS) $(LIBOBJS) $(LIBMIN) $(LIBS)
	- ( cd ../../bin ; ln -fs ../GWW/head/head.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -fv *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

include make.depend
