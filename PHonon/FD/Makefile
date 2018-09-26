# Makefile for school

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3

FDOBJS = \
        stop_pp.o


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a \
         ../../KS_Solvers/libks_solvers.a \
         ../../LAXlib/libqela.a ../../GScratch/libgscratch.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a
PWOBJS = ../../PW/src/libpw.a


TLDEPS= bindir libs pw

all : tldeps fd.x fd_ifc.x fd_ef.x

fd_ifc.x : fd_ifc.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	fd_ifc.o libfd.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

fd.x : fd.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	fd.o libfd.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

fd_ef.x : fd_ef.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	fd_ef.o libfd.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )
tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libfd.a : $(FDOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libfd.a

include make.depend
# DO NOT DELETE
