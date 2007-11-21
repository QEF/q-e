# Makefile for GIPAW

include ../make.sys

GIPAWOBJS = \
gipaw_module.o \
gipaw_main.o \
paw_gipaw.o \
stop_code.o \
apply_p.o \
apply_vel.o \
greenfunction.o \
h_psiq.o \
cg_psi.o \
cgsolve_all.o \
sym_cart_tensor.o \
symmetrize_field.o \
ch_psi_all.o \
test_sum_rule.o \
compute_u_kq.o \
suscept_crystal.o \
j_para.o \
biot_savart.o \
compute_sigma.o \
init_gipaw_1.o \
init_gipaw_2.o \
init_us_2_no_phase.o \
g_tensor_crystal.o \
write_tensor_field.o \
efg.o


MODULES = \
../Modules/atom.o \
../Modules/autopilot.o \
../Modules/basic_algebra_routines.o \
../Modules/berry_phase.o \
../Modules/bfgs_module.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/descriptors.o \
../Modules/dspev_drv.o \
../Modules/electrons_base.o \
../Modules/exc_t.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/input_parameters.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/ions_nose.o \
../Modules/kind.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/mp_base.o \
../Modules/metagga.o \
../Modules/metadyn_base.o \
../Modules/metadyn_io.o \
../Modules/metadyn_vars.o \
../Modules/path_base.o \
../Modules/path_formats.o \
../Modules/path_variables.o \
../Modules/path_opt_routines.o \
../Modules/path_io_routines.o \
../Modules/path_reparametrisation.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/paw_variables.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/radial_grids.o \
../Modules/random_numbers.o \
../Modules/read_cards.o \
../Modules/read_namelists.o \
../Modules/read_upf.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/upf_to_internal.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/vxc_t.o \
../Modules/vxcgc.f90 \
../Modules/wavefunctions.o \
../Modules/wave_base.o \
../Modules/xml_io_base.o \
../Modules/zhpev_drv.o 

PWOBJS = ../PW/libpw.a

TLDEPS=bindir pw mods libs libiotk

all : tldeps gipaw.x

gipaw.x : $(GIPAWOBJS) $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		$(GIPAWOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin; ln -fs ../GIPAW/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f gipaw.x *.o *.F90 *.d *.mod *.i *.L

include make.depend

# DO NOT DELETE

