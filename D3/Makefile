# Makefile for 3rd derivative calculations - D3

include ../make.sys

D3OBJS = \
allocate_d3.o \
bcast_d3_input.o \
ch_psi_all2.o \
close_open.o \
d0rhod2v.o \
d2mxc.o \
d3_exc.o \
d3_init.o \
d3_readin.o \
d3_recover.o \
d3_setup.o \
d3_summary.o \
d3_symdyn.o \
d3_symdynph.o \
d3_valence.o \
d3com.o \
d3dyn_cc.o \
d3ionq.o \
d3matrix.o \
d3toten.o \
d3vrho.o \
davcio_drho2.o \
dpsi_corr.o \
dpsidpsidv.o \
dpsidvdpsi.o \
dqrhod2v.o \
drho_cc.o \
drho_drc.o \
drhod2v.o \
dvdpsi.o \
dvscf.o \
gen_dpdvp.o \
gen_dwfc.o \
incdrhoscf2.o \
openfild3.o \
print_clock_d3.o \
psymd0rho.o \
qstar_d3.o \
read_ef.o \
rotate_and_add_d3.o \
set_d3irr.o \
set_efsh.o \
set_sym_irr.o \
sgama_d3.o \
solve_linter_d3.o \
stop_d3.o \
sym_def1.o \
symd0rho.o \
trntnsc_3.o \
w_1gauss.o \
write_aux.o \
write_d3dyn.o \
write_igk.o \
writed3dyn_5.o

PHOBJS = ../PH/libph.a

MODULES = \
../Modules/atom.o \
../Modules/basic_algebra_routines.o \
../Modules/cell_base.o \
../Modules/clocks.o \
../Modules/check_stop.o \
../Modules/constants.o \
../Modules/control_flags.o \
../Modules/constraints_module.o \
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
../Modules/task_groups.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/upf_to_internal.o \
../Modules/uspp.o \
../Modules/timestep.o \
../Modules/version.o \
../Modules/vxc_t.o \
../Modules/vxcgc.o \
../Modules/upf.o \
../Modules/wavefunctions.o \
../Modules/write_upf_v2.o \
../Modules/xml_io_base.o \
../Modules/zhpev_drv.o

EEMODS=../EE/ee_mod.o ../EE/gcoarse_mod.o

PWOBJS = ../PW/libpw.a

TLDEPS= bindir mods libs pw ph

all : tldeps d3.x

d3.x : $(D3OBJS) $(PHOBJS) $(PWOBJS) $(MODULES) $(EEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o d3.x $(D3OBJS) $(MODULES) $(EEMODS) $(PHOBJS) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../D3/d3.x . )
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f d3.x *.o *~ *.F90 *.d *.i *.mod *.L

include make.depend
