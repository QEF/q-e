#/a Makefile for Modules

include ../make.inc

# location of needed modules
MODFLAGS=$(BASEMOD_FLAGS)

# list of modules

MODULES = \
additional_kpoints.o \
atomic_wfc_mod.o \
autopilot.o \
ions_nose.o \
cell_nose.o \
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
environ_base_module.o \
environment.o \
extffield.o \
fd_gradient.o \
fft_base.o \
fft_rho.o \
fft_wave.o \
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
makov_payne.o \
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
plugin_flags.o \
plugin_arguments.o \
plugin_variables.o \
printout_base.o \
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
copy_order_um_from_xml.o \
setqf.o \
timestep.o\
tsvdw.o\
mbdlib.o\
version.o \
wannier_gw.o\
wannier_new.o \
wavefunctions.o \
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
zvscal.o \
wave_gauge.o

# list of RISM's modules

RISMLIB = \
allocate_fft_3drism.o \
chempot.o \
chempot_lauerism.o \
closure.o \
corrdipole_laue.o \
correctat0_vv.o \
corrgxy0_laue.o \
cryst_to_car_2d.o \
data_structure_3drism.o \
do_1drism.o \
do_3drism.o \
do_lauerism.o \
eqn_1drism.o \
eqn_3drism.o \
eqn_lauedipole.o \
eqn_lauegxy0.o \
eqn_lauelong.o \
eqn_lauerism.o \
eqn_laueshort.o \
eqn_lauevoid.o \
err_rism.o \
guess_3drism.o \
init_1drism.o \
init_3drism.o \
input_1drism.o \
input_3drism.o \
io_rism_xml.o \
lauefft.o \
lauefft_subs.o \
lj_forcefield.o \
lj_solute.o \
molecorr_vv.o \
molebridge_vv.o \
molecule_const.o \
molecule_types.o \
mp_rism.o \
mp_swap_ax_rism.o \
normalize_lauerism.o \
plot_rism.o \
potential_3drism.o \
potential_esm.o \
potential_vv.o \
print_chempot_3drism.o \
print_chempot_lauerism.o \
print_chempot_vv.o \
print_corr_vv.o \
print_solvavg.o \
radfft.o \
read_mol.o \
read_solv.o \
recvec_3drism.o \
rism.o \
rism1d_facade.o \
rism3d_facade.o \
rms_residual.o \
scale_fft_3drism.o \
scale_fft_lauerism.o \
solute.o \
solvation_3drism.o \
solvation_esm.o \
solvation_force.o \
solvation_lauerism.o \
solvation_pbc.o \
solvation_stress.o \
solvavg.o \
solvmol.o \
summary_1drism.o \
summary_3drism.o \
suscept_g0.o \
suscept_laue.o \
suscept_laueint.o \
suscept_vv.o \
write_rism_type.o \
xml_io_rism.o

# list of subroutines and functions (not modules) previously found in flib/clib

OBJS = \
atom_weight.o \
capital.o \
cryst_to_car.o \
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
sort.o \
trimcheck.o \
test_input_file.o \
date_and_tim.o \
volume.o \
wgauss.o \
w0gauss.o \
w1gauss.o \
deviatoric.o \
customize_signals.o \
qmmm_aux.o \
sockets.o

# GPU versions of modules
MODULES += \
  random_numbers_gpu.o

TLDEPS= libfox libla libfft libutil libmbd librxc libupf

all : libqemod.a

libqemod.a: 	$(MODULES) $(OBJS) $(RISMLIB)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend
