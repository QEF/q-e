# Makefile for ZG.x, disca.x, pp_spctrlfn.x, bands_unfold.x, pp_disca.x  

include ../../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../../PW/src \
          $(MOD_FLAG)../../../LR_Modules \
          $(MOD_FLAG)../../../PHonon/PH

PHAUXOBJS = ../../../PHonon/PH/libphaux.a
PHOBJS = ../../../PHonon/PH/libph.a
PWOBJS = ../../../PW/src/libpw.a ../../../KS_Solvers/libks_solvers.a ../../../dft-d3/libdftd3qe.a
LRMODS = ../../../LR_Modules/liblrmod.a
QEMODS = $(BASEMODS)

all : ZG.x disca.x pp_disca.x pp_spctrlfn.x bands_unfold.x epsilon_Gaus.x 

ZG.x : ZG.o $(PHAUXOBJS) $(PHOBJS) $(PWOBJS) $(LRMODS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	ZG.o $(PHAUXOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS) 
	- ( cd ../../../bin ; ln -fs ../EPW/ZG/src/$@ . )

disca.x : disca.o $(PHAUXOBJS) $(PHOBJS) $(PWOBJS) $(LRMODS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	disca.o $(PHAUXOBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS) 
	- ( cd ../../../bin ; ln -fs ../EPW/ZG/src/$@ . )

pp_disca.x : pp_disca.o $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	pp_disca.o $(PHAUXOBJS) $(PHOBJS) $(PWOBJS) $(QEMODS) $(QELIBS) 
	- ( cd ../../../bin ; ln -fs ../EPW/ZG/src/$@ . )

bands_unfold.x : bands_unfold.o $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	bands_unfold.o ../../../PP/src/libpp.a $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../../bin ; ln -fs ../EPW/ZG/src/$@ . )

pp_spctrlfn.x : pp_spctrlfn.o $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	pp_spctrlfn.o $(PHAUXOBJS) $(PHOBJS) $(PWOBJS) $(QEMODS) $(QELIBS) 
	- ( cd ../../../bin ; ln -fs ../EPW/ZG/src/$@ . )

epsilon_Gaus.x : epsilon_Gaus.o $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	epsilon_Gaus.o ../../../PP/src/libpp.a $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../../bin ; ln -fs ../EPW/ZG/src/$@ . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *.mod
	- /bin/rm -f ../../../bin/ZG.x
	- /bin/rm -f ../../../bin/disca.x
	- /bin/rm -f ../../../bin/pp_disca.x
	- /bin/rm -f ../../../bin/bands_unfold.x
	- /bin/rm -f ../../../bin/pp_spctrlfn.x
	- /bin/rm -f ../../../bin/epsilon_Gaus.x
	- /bin/rm -f local/*.x local/*.o

include make.depend
