# Makefile for tools

include ../../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../src \
          $(MOD_FLAG)../../dft-d3/
QEMODS = ../../dft-d3/libdftd3qe.a $(BASEMODS)
PWOBJS = ../src/libpw.a ../../KS_Solvers/libks_solvers.a

# see comments in file tldeps for explanations
TLDEPS= pwlibs
include ../../install/tldeps

all : ev.x kpoints.x pwi2xsf.x ibrav2cell.x cell2ibrav.x scan_ibrav.x rism1d.x

ibrav2cell.x : ibrav2cell.o $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		ibrav2cell.o $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )
cell2ibrav.x : cell2ibrav.o $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		cell2ibrav.o $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )
scan_ibrav.x : scan_ibrav.o $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		scan_ibrav.o $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

ev.x : ev.o $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		ev.o $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

kpoints.x : kpoints.o $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		kpoints.o $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

pwi2xsf.x : pwi2xsf.o $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		pwi2xsf.o $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

rism1d.x : rism1d.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		rism1d.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.mod *.d *.i *.L


include make.depend


