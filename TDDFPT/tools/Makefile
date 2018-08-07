# Makefile for TDDFPT tools

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../KS_Solvers/CG \
          $(MOD_FLAG)../../KS_Solvers/PPCG \
          $(MOD_FLAG)../../KS_Solvers/Davidson \
          $(MOD_FLAG)../../LR_Modules

QEMODS = ../../Modules/libqemod.a \
         ../../KS_Solvers/CG/libcg.a ../../KS_Solvers/PPCG/libppcg.a ../../KS_Solvers/Davidson/libdavid.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a
PWOBJS = ../../PW/src/libpw.a

TLDEPS=bindir libs pw

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

