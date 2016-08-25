# Makefile for school

include ../../make.inc

# location of include files
#IFLAGS=-I../../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG).

FDOBJS = \
        stop_pp.o 


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a


TLDEPS= bindir libs pw

all : tldeps fd.x fd_ifc.x fd_ef.x

fd_ifc.x : fd_ifc.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
	fd_ifc.o libfd.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

fd.x : fd.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
	fd.o libfd.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) 
	- ( cd ../../bin ; ln -fs ../PHonon/FD/$@ . )

fd_ef.x : fd_ef.o libfd.a  $(FDOBJS) $(PWOBJS) $(QEMODS)  
	$(LD) $(LDFLAGS) -o $@ \
	fd_ef.o libfd.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
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

