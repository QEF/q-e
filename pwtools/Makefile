# Makefile for tools

include ../make.sys
include ../make.rules

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
../Modules/io_global.o \
../Modules/io_files.o \
../Modules/fft_scalar.o \
../Modules/kind.o \
../Modules/mp.o \
../Modules/mp_global.o \
../Modules/parallel_include.o \
../Modules/shmem_include.o \
../Modules/parameters.o \
../Modules/supercell.o \
../Modules/splinelib.o

all : band_plot.x dist.x dynmat.x ev.x fqha.x kpoints.x matdyn.x \
      path_int.x pwi2xsf.x q2r.x

# dependencies
band_plot.x : band_plot.o
fqha.x : fqha.o

dist.x : dist.o
ev.x : ev.o
kpoints.x : kpoints.o
lambda.x : lambda.o
path_int.x : path_int.o
pwi2xsf.x : pwi2xsf.o

dynmat.x : dynmat.o rigid.o
matdyn.x : matdyn.o rigid.o
q2r.x : q2r.o rigid.o

# compilation rules
band_plot.x fqha.x :
	$(LD) -o $@ $<  $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dist.x ev.x kpoints.x lambda.x path_int.x pwi2xsf.x :
	$(LD) -o $@ $< $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

dynmat.x matdyn.x q2r.x :
	$(LD) -o $@ $< rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../pwtools/$@ . )

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include .dependencies
