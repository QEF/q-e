# Makefile for TDDFPT tools

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../LR_Modules

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a \
	 ../../upflib/libupf.a    ../../KS_Solvers/libks_solvers.a \
         ../../LAXlib/libqela.a ../../GScratch/libgscratch.a ../../UtilXlib/libutil.a
PWOBJS = ../../PW/src/libpw.a

TLDEPS=pwlibs

all : tldeps turbo_spectrum.x

turbo_spectrum.x : tddfpt_calculate_spectrum.o $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	tddfpt_calculate_spectrum.o $(QEMODS) $(PWOBJS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/tools/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../tools/$@ . ); fi

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	rm -f ../bin/turbo_spectrum.x *.x *.o *.mod *_tmp.f90
