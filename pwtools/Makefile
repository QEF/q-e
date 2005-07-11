# Makefile for tools

include ../make.sys

PWOBJS = \
../PW/cft_3.o \
../PW/checksym.o \
../PW/coset.o \
../PW/cryst_to_car.o \
../PW/cubicsym.o \
../PW/eqvect.o \
../PW/error_handler.o \
../PW/hexsym.o \
../PW/irrek.o \
../PW/kpoint_grid.o \
../PW/mode_group.o \
../PW/multable.o \
../PW/noncol.o \
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
../Modules/fft_scalar.o \
../Modules/kind.o \
../Modules/mp.o \
../Modules/mp_global.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/path_formats.o \
../Modules/splinelib.o \
../Modules/timestep.o

all : band_plot.x dist.x dynmat.x ev.x fqha.x kpoints.x matdyn.x \
      path_int.x pwi2xsf.x q2r.x

band_plot.x : band_plot.o
	$(LD) -o $@ band_plot.o $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dist.x : dist.o
	$(LD) -o $@ dist.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dynmat.x : dynmat.o rigid.o
	$(LD) -o $@ dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

ev.x : ev.o
	$(LD) -o $@ ev.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

fqha.x : fqha.o
	$(LD) -o $@ fqha.o $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

kpoints.x : kpoints.o
	$(LD) -o $@ kpoints.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

lambda.x : lambda.o
	$(LD) -o $@ lambda.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

matdyn.x : matdyn.o rigid.o
	$(LD) -o $@ matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

path_int.x : path_int.o
	$(LD) -o $@ path_int.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

pwi2xsf.x : pwi2xsf.o
	$(LD) -o $@ pwi2xsf.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

q2r.x : q2r.o rigid.o
	$(LD) -o $@ q2r.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
