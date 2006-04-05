# Makefile for tools

include ../make.sys

PWOBJS = \
../PW/checksym.o \
../PW/coset.o \
../PW/cryst_to_car.o \
../PW/cubicsym.o \
../PW/eqvect.o \
../PW/error_handler.o \
../PW/hexsym.o \
../PW/irrek.o \
../PW/inverse_s.o \
../PW/kpoint_grid.o \
../PW/mode_group.o \
../PW/multable.o \
../PW/noncol.o \
../PW/pwcom.o \
../PW/sgam_at.o \
../PW/sgam_at_mag.o \
../PW/sgam_ph.o \
../PW/sgama.o \
../PW/smallg_q.o \
../PW/trnvecc.o \
../PW/w0gauss.o \
../PW/wsweight.o

MODULES = \
../Modules/basic_algebra_routines.o \
../Modules/cell_base.o \
../Modules/constants.o \
../Modules/control_flags.o \
../Modules/io_global.o \
../Modules/io_files.o \
../Modules/ions_base.o \
../Modules/fft_types.o \
../Modules/fft_scalar.o \
../Modules/kind.o \
../Modules/mp.o \
../Modules/mp_global.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/path_formats.o \
../Modules/printout_base.o \
../Modules/random_numbers.o \
../Modules/recvec.o \
../Modules/splinelib.o \
../Modules/timestep.o \
../Modules/uspp.o

TLDEPS= bindir mods libs pw

all : tldeps band_plot.x dist.x dynmat.x ev.x fqha.x kpoints.x matdyn.x \
      path_int.x pwi2xsf.x q2r.x bands_FS.x kvecs_FS.x 

band_plot.x : band_plot.o $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ band_plot.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dist.x : dist.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		dist.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dynmat.x : dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

ev.x : ev.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		ev.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

fqha.x : fqha.o $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ fqha.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

kpoints.x : kpoints.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		kpoints.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

lambda.x : lambda.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		lambda.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

matdyn.x : matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

metadyn_pp.x : metadyn_pp.o $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ metadyn_pp.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

path_int.x : path_int.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		path_int.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

pwi2xsf.x : pwi2xsf.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pwi2xsf.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

q2r.x : q2r.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		q2r.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )
# Isaev
bands_FS.x : bands_FS.o
	$(MPIF90) $(LDFLAGS) -o $@ bands_FS.o $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

kvecs_FS.x : kvecs_FS.o
	$(MPIF90) $(LDFLAGS) -o $@ kvecs_FS.o $(LIBS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
