# Makefile for phonon (PH)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)
TLDEPS= phlibs
PHOBJS = phq_readin2.o do_phonon2.o supp.o ../../PHonon/PH/libph.a

all : ph.x

ph.x : phonon2.o $(PHOBJS) $(PWOBJS) $(LRMODS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ phonon2.o $(PHOBJS) \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
#	- ( cd ../../bin ; ln -fs ../PHonon/PH/ph.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

depend:
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)

phonon2.o: supp.o
do_phonon2.o: supp.o
phq_readin2.o: supp.o
include ../../PHonon/PH/make.depend
