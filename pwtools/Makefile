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
../Modules/formats.o \
../Modules/kind.o \
../Modules/mp.o \
../Modules/mp_global.o \
../Modules/parallel_include.o \
../Modules/shmem_include.o \
../Modules/parameters.o \
../Modules/supercell.o \
../Modules/splinelib.o

all : band_plot.x q2r.x matdyn.x dynmat.x fqha.x dist.x ev.x \
      kpoints.x pwi2xsf.x path_int.x

band_plot.x : band_plot.o
	$(LD) -o band_plot.x band_plot.o  $(LDFLAGS)

dist.x : dist.o
	$(LD) -o dist.x dist.o $(PWOBJS) $(MODULES) $(LDFLAGS)

q2r.x : q2r.o rigid.o
	$(LD) -o q2r.x q2r.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)

matdyn.x : matdyn.o rigid.o
	$(LD) -o matdyn.x matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)

dynmat.x : dynmat.o rigid.o
	$(LD) -o dynmat.x dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)

fqha.x : fqha.o
	$(LD) -o fqha.x fqha.o  $(LDFLAGS)

lambda.x : lambda.o
	$(LD) -o lambda.x lambda.o $(PWOBJS) $(MODULES) $(LDFLAGS)

ev.x : ev.o
	$(LD) -o ev.x ev.o $(PWOBJS) $(MODULES) $(LDFLAGS)

kpoints.x : kpoints.o
	$(LD) -o kpoints.x kpoints.o $(PWOBJS) $(MODULES) $(LDFLAGS)

pwi2xsf.x : pwi2xsf.o
	$(LD) -o pwi2xsf.x pwi2xsf.o $(PWOBJS) $(MODULES) $(LDFLAGS)

path_int.x : path_int.o
	$(LD) -o path_int.x path_int.o $(PWOBJS) $(MODULES) $(LDFLAGS)

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include .dependencies
