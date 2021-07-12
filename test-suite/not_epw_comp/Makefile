# Makefile for phonon (PH)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a  ../../XClib/xc_lib.a ../../upflib/libupf.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a
TLDEPS= phlibs
PHOBJS = phq_readin2.o do_phonon2.o supp.o ../../PHonon/PH/libph.a

all : ph.x

ph.x : phonon2.o $(PHOBJS) $(PWOBJS) $(LRMODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ phonon2.o $(PHOBJS) \
		 $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
#	- ( cd ../../bin ; ln -fs ../PHonon/PH/ph.x . )

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

depend:
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)

include make.depend
# DO NOT DELETE
