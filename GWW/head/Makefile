# Makefile for head (head.x)
# 
# Author G. Stenuit
#
include ../../make.sys

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PH $(MOD_FLAG)../../PW \
	  $(MOD_FLAG)../pw4gww $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../flib/flib.a \
         ../../clib/clib.a ../../flib/ptools.a

HEADOBJS = \
bcast_ph_input.o \
bcast_ph_input1.o \
close_phq.o \
openfilq.o \
phq_init.o \
phq_readin.o \
solve_head.o 

QEMODS = ../../Modules/libqemod.a
LIBPWPH = ../pw4gww/libpw4gww.a ../../PH/libph.a ../../PW/libpw.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS= bindir mods libs pw ph 

all : tldeps head.x  

head.x : head.o $(LIBOBJS) $(HEADOBJS) $(LIBPWPH) $(LIBMIN)
	$(LD) $(LDFLAGS) -o head.x head.o \
		 $(HEADOBJS) $(LIBPWPH) $(QEMODS) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../GWW/head/head.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

include make.depend
