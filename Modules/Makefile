#/a Makefile for Modules

include ../make.inc

# location of needed modules
MODFLAGS=$(BASEMOD_FLAGS) \
         $(MOD_FLAG)../ELPA/src

# list of modules

MODULES = \
additional_kpoints.o \
autopilot.o \
basic_algebra_routines.o \
becmod.o \
bfgs_module.o \
bspline.o \
bz_form.o \
cell_base.o  \
check_stop.o  \
command_line_options.o \
compute_dipole.o \
constants.o \
constraints_module.o \
control_flags.o \
coulomb_vcut.o \
dist.o \
electrons_base.o \
environment.o \
fd_gradient.o \
fft_base.o \
fft_rho.o \
fsockets.o \
funct.o \
generate_function.o \
gradutils.o \
gvecw.o \
input_parameters.o \
invmat.o \
io_files.o \
io_global.o  \
ions_base.o \
kind.o \
lmdif.o \
mdiis.o \
mm_dispersion.o \
mp_bands.o \
mp_exx.o \
mp_global.o \
mp_images.o \
mp_pools.o \
mp_wave.o \
mp_world.o \
noncol.o \
open_close_input_file.o \
parameters.o \
parser.o \
paw_variables.o \
plugin_flags.o \
plugin_arguments.o \
plugin_variables.o \
pw_dot.o \
qmmm.o \
random_numbers.o \
read_cards.o \
read_input.o \
read_namelists.o \
read_pseudo.o \
recvec.o \
recvec_subs.o \
run_info.o \
space_group.o \
set_para_diag.o \
set_signal.o \
set_vdw_corr.o \
setqf.o \
timestep.o\
tsvdw.o\
mbdlib.o\
version.o \
wannier_gw.o\
wannier_new.o \
wavefunctions.o \
wrappers.o\
ws_base.o \
xc_vdW_DF.o \
xc_rVV10.o \
io_base.o \
qes_types_module.o \
qes_libs_module.o  \
qes_write_module.o \
qes_read_module.o \
qes_reset_module.o \
qes_init_module.o \
qes_read_module.o \
qes_bcast_module.o \
qexsd.o   \
qexsd_copy.o   \
qexsd_init.o   \
qexsd_input.o \
hdf5_qe.o\
qeh5_module.o\
fox_init_module.o \
xsf.o \
wyckoff.o \
wypos.o \
zvscal.o 

# list of subroutines and functions (not modules) previously found in flib/

OBJS = \
atom_weight.o \
capital.o \
cryst_to_car.o \
erf.o \
expint.o \
generate_k_along_lines.o \
has_xml.o \
inpfile.o \
int_to_char.o \
latgen.o \
linpack.o \
matches.o \
plot_io.o \
radial_gradients.o \
rgen.o \
recips.o \
remove_tot_torque.o \
set_hubbard_l.o \
set_hubbard_n.o \
sort.o \
sph_bes.o \
sph_dbes.o \
trimcheck.o \
test_input_file.o \
date_and_tim.o \
volume.o \
dylmr2.o \
wgauss.o \
w0gauss.o \
w1gauss.o \
deviatoric.o

# GPU versions of modules
MODULES += \
  wavefunctions_gpu.o \
  recvec_gpu.o \
  becmod_gpu.o \
  becmod_subs_gpu.o \
  cuda_subroutines.o \
  dylmr2_gpu.o \
  random_numbers_gpu.o \
  sph_bes_gpu.o

TLDEPS= libfox libutil libla libfft librxc

all : libqemod.a

## The following is needed only for lapack compiled from sources

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

libqemod.a: 	$(MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend
