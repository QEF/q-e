# Makefile for ZG.x

include ../../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../../PW/src \
          $(MOD_FLAG)../../../LR_Modules \
          $(MOD_FLAG)../../../PHonon/PH

PHAUXOBJS = ../../../PHonon/PH/libphaux.a
PHOBJS = ../../../PHonon/PH/libph.a
PWOBJS = ../../../PW/src/libpw.a
LRMODS = ../../../LR_Modules/liblrmod.a
QEMODS = ../../../Modules/libqemod.a ../../../KS_Solvers/libks_solvers.a \
         ../../../upflib/libupf.a \
         ../../../FFTXlib/libqefft.a ../../../LAXlib/libqela.a ../../../UtilXlib/libutil.a ../../../dft-d3/libdftd3qe.a
LIBOBJS = ../../../clib/clib.a

ZG.x : ZG.o $(PHAUXOBJS) $(PHOBJS) $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	ZG.o $(PHAUXOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) 
	- ( cd ../../../bin ; ln -fs ../EPW/ZG_displacement/src/$@ . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *.mod
	- /bin/rm -f ../../../bin/ZG.x

include make.depend
