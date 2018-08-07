# Makefile for tools

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../KS_Solvers/Davidson \
          $(MOD_FLAG)../../KS_Solvers/CG \
          $(MOD_FLAG)../../KS_Solvers/PPCG \
          $(MOD_FLAG)../src
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a \
         ../../KS_Solvers/CG/libcg.a ../../KS_Solvers/PPCG/libppcg.a ../../KS_Solvers/Davidson/libdavid.a \
         ../../LAXlib/libqela.a ../../UtilXlib/libutil.a
PWOBJS = ../src/libpw.a

TLDEPS= bindir mods libs pw libfft libla libutil

LIBOBJS        = ../../clib/clib.a ../../iotk/src/libiotk.a

all : tldeps ev.x kpoints.x pwi2xsf.x ibrav2cell.x

ibrav2cell.x : ibrav2cell.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		ibrav2cell.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )
	- ( cd ../../bin ; ln -fs ../PW/tools/cell2ibrav.py cell2ibrav.x )

ev.x : ev.o ev_xml.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		ev.o ev_xml.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

kpoints.x : kpoints.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		kpoints.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

pwi2xsf.x : pwi2xsf.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pwi2xsf.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

tldeps:
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *_tmp.f90 *.mod *.d *.i *.L

include make.depend
