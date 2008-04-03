# Makefile for Gamma-only code

include ../make.sys

PWOBJS = ../PW/libpw.a

MODULES = \
../Modules/atom.o \
../Modules/bfgs_module.o \
../Modules/basic_algebra_routines.o \
../Modules/berry_phase.o \
../Modules/cell_base.o \
../Modules/clocks.o \
../Modules/check_stop.o \
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
../Modules/input_parameters.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/kind.o \
../Modules/metadyn_base.o \
../Modules/metadyn_io.o \
../Modules/metadyn_vars.o \
../Modules/metagga.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/mp_base.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/path_base.o \
../Modules/path_formats.o \
../Modules/path_variables.o \
../Modules/path_opt_routines.o \
../Modules/path_io_routines.o \
../Modules/path_reparametrisation.o \
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
../Modules/task_groups.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/timestep.o \
../Modules/upf_to_internal.o \
../Modules/upf.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/vxc_t.o \
../Modules/vxcgc.o \
../Modules/wavefunctions.o \
../Modules/write_upf_v2.o \
../Modules/xml_io_base.o \
../Modules/zhpev_drv.o

CGOBJS = \
a_h.o \
cg_readin.o \
cg_setup.o \
cg_setupdgc.o \
cgcom.o \
cgsolve.o \
d2ion.o \
dgcxc.o \
dgcxc_spin.o \
dgradcorr.o \
dielec.o \
drhodv.o \
dvb_cc.o \
dvpsi_e.o \
dvpsi_kb.o \
dyndiar.o \
dynmat_init.o \
dynmatcc.o \
find_equiv_sites.o \
generate_dynamical_matrix.o \
generate_effective_charges.o \
h_h.o \
macro.o \
pw_dot.o \
phcg.o \
rhod2vkb.o \
solve_e.o \
solve_ph.o \
writedyn.o \
zvscal.o

TLDEPS= bindir mods libs pw

all : tldeps phcg.x

phcg.x : $(PWOBJS) $(CGOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o phcg.x \
		$(CGOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../Gamma/phcg.x . )
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *.F90 *.mod *.L

include make.depend
