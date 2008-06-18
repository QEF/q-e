# Makefile for TF-VdW calculation.

include ../make.sys


VdWOBJS = \
allocate_vdw.o \
cegterg.o \
check_v_eff.o \
ch_psi_all.o \
dv_of_drho.o \
dvpsi_e.o \
eff_pot.o \
h_psiq.o \
incdrhoscf.o \
pbcg_psi.o \
polariz.o \
print_clock_vdw.o \
solve_e.o \
start_vdw.o \
stop_vdw.o \
vdwcom.o

PPOBJS = \
../PP/openfil_pp.o

MODULES = \
../Modules/atom.o \
../Modules/basic_algebra_routines.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/descriptors.o \
../Modules/dspev_drv.o \
../Modules/electrons_base.o \
../Modules/error_handler.o \
../Modules/exc_t.o \
../Modules/fft_base.o \
../Modules/fft_parallel.o \
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
../Modules/mp_base.o \
../Modules/metagga.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/paw_variables.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/radial_grids.o \
../Modules/random_numbers.o \
../Modules/read_ncpp.o \
../Modules/read_upf_v1.o \
../Modules/read_upf_v2.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/stick_base.o \
../Modules/splinelib.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/upf.o \
../Modules/upf_to_internal.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/vxcgc.o \
../Modules/vxc_t.o \
../Modules/wavefunctions.o \
../Modules/write_upf_v2.o \
../Modules/xml_io_base.o \
../Modules/zhpev_drv.o

EEMODS=../EE/ee_mod.o ../EE/gcoarse_mod.o

PWOBJS = ../PW/libpw.a

PHOBJS = ../PH/libph.a

TLDEPS= bindir mods libs pw ph pp

all : tldeps vdw.x
 
vdw.x : vdw.o $(VdWOBJS) $(PPOBJS) $(PWOBJS) $(PHOBJS) $(MODULES) $(EEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		vdw.o $(VdWOBJS) $(MODULES) $(EEMODS) $(PPOBJS) $(PHOBJS) $(PWOBJS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../VdW/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L

include make.depend
