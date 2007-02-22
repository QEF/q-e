# Makefile for postprocessing (PP)

include ../make.sys


PPOBJS = \
addusdens1d.o \
add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
atomic_wfc_nc_proj.o \
cft.o \
cgracsc.o \
compute_ppsi.o \
compute_sigma_avg.o \
cube.o \
cubicspinsym.o \
dosg.o \
do_initial_state.o \
do_shift_ew.o \
d_matrix_nc.o \
d_matrix_so.o \
elf.o \
ggen1d.o \
hexspinsym.o \
local_dos.o \
local_dos_mag.o \
local_dos1d.o \
openfil_pp.o\
plot_io.o \
plot_whole_cell.o \
punch_plot.o \
smallgk.o \
start_postproc.o \
stm.o \
stop_pp.o \
sym_band.o \
work_function.o \
write_p_avg.o \
xsf.o 

MODULES = \
../Modules/atom.o \
../Modules/basic_algebra_routines.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/electrons_base.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/input_parameters.o \
../Modules/kind.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/metagga.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/random_numbers.o \
../Modules/read_upf.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/stick_base.o \
../Modules/splinelib.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o \
../Modules/wannier.o \
../Modules/xml_io_base.o

PWOBJS = ../PW/libpw.a

TLDEPS= bindir mods libs pw

all : tldeps average.x bands.x dos.x plotband.x plotproj.x \
      plotrho.x wfdd.x \
      pmw.x pp.x projwfc.x pw2casino.x pw2wannier90.x pw_export.x \
      voronoy.x initial_state.x dipole.x sumpdos.x plan_avg.x pw2gw.x

average.x : average.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		average.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

bands.x : bands.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		bands.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

dipole.x : dipole.o $(PPOBJS) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		dipole.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

dos.x : dos.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		dos.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plan_avg.x : plan_avg.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plan_avg.o $(PPOBJS) $(MODULES) $(PWOBJS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotband.x : plotband.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plotband.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotproj.x : plotproj.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plotproj.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotrho.x : plotrho.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plotrho.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pmw.x : poormanwannier.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		poormanwannier.o $(PPOBJS) $(MODULES) $(PWOBJS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

xctest.x : xctest.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		xctest.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pp.x : postproc.o chdens.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		postproc.o chdens.o $(PPOBJS) $(MODULES) $(PWOBJS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

projwfc.x : projwfc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		projwfc.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2casino.x : pw2casino.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2casino.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2wannier90.x : pw2wannier90.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2wannier90.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

voronoy.x : voronoy.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		voronoy.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

initial_state.x : initial_state.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		initial_state.o $(PPOBJS) $(MODULES) $(PWOBJS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2gw.x : pw2gw.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2gw.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw_export.x : pw_export.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw_export.o $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

wfdd.x : wfdd.o
	$(MPIF90) $(LDFLAGS) -o $@ wfdd.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

sumpdos.x : sumpdos.o 
	$(MPIF90) $(LDFLAGS) -o $@ sumpdos.o
	- ( cd ../bin ; ln -fs ../PP/$@ . )

epsilon.x : epsilon.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ epsilon.o $(PPOBJS) $(MODULES) $(PWOBJS) \
	$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L

include make.depend
