# Makefile for postprocessing (PP)

include ../make.sys

PPOBJS = \
../PP/start_postproc.o \
../PP/plot_io.o \
../PP/stop_pp.o

PWCONDOBJS = \
allocate_cond.o \
bessj.o \
compbs.o \
compbs_2.o \
condcom.o \
condmain.o \
cond_out.o \
do_cond.o \
eigenchnl.o \
form_zk.o \
four.o \
free_mem.o \
gep_x.o \
gramsh.o \
hev_ab.o \
init_cond.o \
init_orbitals.o \
init_gper.o \
integrals.o \
jbloch.o \
kbloch.o \
local.o \
local_set.o \
poten.o \
print_clock_pwcond.o \
rotproc.o \
save_cond.o \
scatter_forw.o \
summary_band.o \
summary_tran.o \
sunitary.o \
transmit.o \
write_states.o

MODULES = \
../Modules/atom.o \
../Modules/basic_algebra_routines.o \
../Modules/berry_phase.o \
../Modules/cell_base.o \
../Modules/clocks.o \
../Modules/check_stop.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/descriptors.o \
../Modules/electrons_base.o \
../Modules/exc_t.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/grid_paw_variables.o \
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
../Modules/radial_grids.o \
../Modules/read_paw.o \
../Modules/read_upf.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/task_groups.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/timestep.o \
../Modules/upf_to_internal.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/vxc_t.o \
../Modules/vxcgc.o \
../Modules/wavefunctions.o \
../Modules/xml_io_base.o


PWOBJS = ../PW/libpw.a

TLDEPS= bindir mods libs pw pp

all : tldeps pwcond.x

pwcond.x : $(PWCONDOBJS) $(PWOBJS) $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o pwcond.x \
		$(PWCONDOBJS) $(PPOBJS) $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PWCOND/pwcond.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :
 
clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.a *.L

include make.depend
