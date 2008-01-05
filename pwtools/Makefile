# Makefile for tools

include ../make.sys

PWOBJS = ../PW/libpw.a

MODULES = \
../Modules/basic_algebra_routines.o \
../Modules/cell_base.o \
../Modules/constants.o \
../Modules/control_flags.o \
../Modules/error_handler.o \
../Modules/io_global.o \
../Modules/io_files.o \
../Modules/ions_base.o \
../Modules/fft_types.o \
../Modules/fft_scalar.o \
../Modules/kind.o \
../Modules/mp.o \
../Modules/mp_base.o \
../Modules/mp_global.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/path_formats.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/random_numbers.o \
../Modules/radial_grids.o \
../Modules/recvec.o \
../Modules/splinelib.o \
../Modules/timestep.o \
../Modules/uspp.o

TLDEPS= bindir mods libs pw

all : tldeps band_plot.x dist.x ev.x fqha.x kpoints.x \
      path_int.x pwi2xsf.x bands_FS.x kvecs_FS.x metadyn_pp.x

band_plot.x : band_plot.o $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ band_plot.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dist.x : dist.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		dist.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

ev.x : ev.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		ev.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

fqha.x : fqha.o $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fqha.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

kpoints.x : kpoints.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		kpoints.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

lambda.x : lambda.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		lambda.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

metadyn_pp.x : metadyn_pp.o $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ metadyn_pp.o $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

path_int.x : path_int.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		path_int.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

pwi2xsf.x : pwi2xsf.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pwi2xsf.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

# Isaev
bands_FS.x : bands_FS.o
	$(LD) $(LDFLAGS) -o $@ bands_FS.o $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

kvecs_FS.x : kvecs_FS.o
	$(LD) $(LDFLAGS) -o $@ kvecs_FS.o $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend
