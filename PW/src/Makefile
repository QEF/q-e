# Makefile for PW

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib

PWOBJS = \
pwscf.o 

PWLIBS = \
a2fmod.o \
add_bfield.o \
add_efield.o \
add_vuspsi.o \
add_monofield.o \
add_paw_to_deeq.o \
add_vhub_to_deeq.o \
addusdens.o \
addusforce.o \
addusstress.o \
allocate_fft.o \
allocate_locpot.o \
allocate_nlpot.o \
allocate_wfc.o \
atomic_rho.o \
atomic_wfc.o \
atomic_wfc_mod.o \
average_pp.o \
acfdt_in_pw.o \
newd.o \
bp_mod.o \
bp_c_phase.o \
bp_calc_btq.o \
bp_qvan3.o \
bp_strings.o \
buffers.o \
c_bands.o \
c_phase_field.o \
orbm_kubo.o \
ccgdiagg.o \
cdiagh.o \
cdiaghg.o \
cegterg.o \
clean_pw.o \
close_files.o \
compute_becsum.o \
compute_deff.o \
compute_dip.o \
compute_rho.o \
compute_qdipol.o \
compute_qdipol_so.o \
compute_ux.o \
coset.o \
d_matrix.o \
data_structure.o \
data_structure_custom.o \
deriv_drhoc.o \
divide_class.o \
divide_class_so.o \
realus.o \
divide.o \
divide_et_impera.o \
dqvan2.o \
drhoc.o \
run_driver.o \
dvloc_of_g.o \
dynamics_module.o \
efermig.o \
efermit.o \
electrons.o \
eqvect.o \
esm.o \
esm_cft.o \
ewald.o \
ewald_dipole.o \
exx.o \
fcp.o \
find_group.o \
forces_bp_efield.o \
force_cc.o \
force_corr.o \
force_ew.o \
force_hub.o \
force_lc.o \
force_us.o \
forces.o \
g_psi.o \
g_psi_mod.o \
gen_at_dj.o \
gen_at_dy.o \
gen_us_dj.o \
gen_us_dy.o \
get_locals.o \
gk_sort.o \
gradcorr.o \
gweights.o \
g2_kin.o \
h_epsi_her_apply.o \
h_epsi_her_set.o \
h_1psi.o \
h_psi.o \
h_psi_meta.o \
hinit0.o \
hinit1.o \
init_ns.o \
init_q_aeps.o \
init_run.o \
init_us_0.o \
init_us_b0.o \
init_us_1.o \
init_us_2.o \
init_at_1.o \
init_vloc.o \
input.o \
interpolate.o \
io_rho_xml.o \
irrek.o \
iweights.o \
start_k.o \
kpoint_grid.o \
lchk_tauxk.o \
ldaU.o \
make_pointlists.o \
makov_payne.o \
manypw.o \
martyna_tuckerman.o \
memory_report.o \
mix_rho.o \
move_ions.o \
multable.o \
n_plane_waves.o \
new_ns.o \
new_occ.o \
ns_adj.o \
non_scf.o \
offset_atom_wfc.o \
openfil.o \
orthoatwfc.o \
output_tau.o \
para.o \
paw_exx.o \
paw_init.o \
paw_onecenter.o \
paw_symmetry.o \
plugin_print_energies.o \
plugin_scf_energy.o \
plugin_scf_potential.o \
plugin_init_ions.o \
plugin_init_cell.o \
plugin_init_potential.o \
plugin_initbase.o \
plugin_clean.o \
plugin_clock.o \
plugin_summary.o \
plugin_initialization.o \
plugin_ext_forces.o \
plugin_int_forces.o \
plugin_read_input.o \
plus_u_full.o \
potinit.o \
print_clock_pw.o \
print_ks_energies.o \
punch.o \
pw_restart.o \
pw_restart_new.o \
add_qexsd_step.o \
pw_init_qexsd_input.o \
pwcom.o \
pw2blip.o \
pw2casino.o \
pw2casino_write.o \
qvan2.o \
rcgdiagg.o \
rdiagh.o \
rdiaghg.o \
read_conf_from_file.o \
read_file.o \
read_file_new.o \
regterg.o \
remove_atomic_rho.o \
report_mag.o \
restart_in_electrons.o \
rho2zeta.o \
rotate_wfc.o \
rotate_wfc_k.o \
rotate_wfc_gamma.o \
ruotaijk.o \
run_pwscf.o \
s_1psi.o \
s_psi.o \
save_in_cbands.o \
save_in_electrons.o \
scale_h.o \
scf_mod.o \
set_kplusq.o \
set_kup_and_kdw.o \
set_rhoc.o \
set_vrs.o \
setlocal.o \
setup.o \
spinor.o \
sph_ind.o \
stop_run.o \
stres_cc.o \
stres_ewa.o \
stres_gradcorr.o \
stres_har.o \
stres_hub.o \
stres_knl.o \
stres_loc.o \
stres_us.o \
stres_nonloc_dft.o \
stress.o \
struct_fact.o \
sum_band.o \
sumkg.o \
sumkt.o \
summary.o \
symme.o \
symm_base.o \
symmetrize_at.o \
tabd.o \
tetra.o \
transform_becsum_so.o \
transform_becsum_nc.o \
transform_qq_so.o \
trnvecc.o \
update_pot.o \
us_exx.o \
usnldiag.o \
v_of_rho.o \
vcsmd.o \
vcsubs.o \
vhpsi.o \
vloc_of_g.o \
vloc_psi.o \
xdm_dispersion.o \
wfcinit.o \
write_ns.o \
wsweight.o \
weights.o \
ortho_wfc.o \
wannier_proj.o \
wannier_init.o \
wannier_check.o \
wannier_clean.o \
wannier_occ.o \
wannier_enrg.o \
wyckoff.o       

QEMODS=../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS=bindir mods libs liblapack libblas

all : tldeps pw.x generate_vdW_kernel_table.x generate_rVV10_kernel_table.x

pw.x : $(PWOBJS) libpw.a $(LIBOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(PWOBJS) libpw.a $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin; ln -fs ../PW/src/$@ . ; \
               ln -fs ../PW/src/$@ dist.x ; ln -fs ../PW/src/$@ manypw.x ; )
# dist.x and manypw.x link the same executable pw.x

libpw.a : $(PWLIBS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

generate_vdW_kernel_table.x : libpw.a generate_vdW_kernel_table.o
	$(LD) $(LDFLAGS) -o $@ \
		generate_vdW_kernel_table.o $(QEMODS) libpw.a $(LIBOBJS) $(LIBS)
	- ( cd ../../bin; ln -fs ../PW/src/$@ . )

generate_rVV10_kernel_table.x : libpw.a generate_rVV10_kernel_table.o
	$(LD) $(LDFLAGS) -o $@ \
		generate_rVV10_kernel_table.o $(QEMODS) libpw.a $(LIBOBJS) $(LIBS)
	- ( cd ../../bin; ln -fs ../PW/src/$@ . )


include make.depend

# DO NOT DELETE
