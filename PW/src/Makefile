# Makefile for PW

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../dft-d3/

PWOBJS = \
pwscf.o 

PWLIBS = \
a2fmod.o \
add_bfield.o \
add_efield.o \
add_dmft_occ.o \
add_vuspsi.o \
add_gatefield.o \
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
beef.o \
bp_mod.o \
bp_c_phase.o \
bp_strings.o \
buffers.o \
c_bands.o \
c_phase_field.o \
orbm_kubo.o \
cdiagh.o \
clean_pw.o \
close_files.o \
commutator_Hx_psi.o \
commutator_Vhubx_psi.o \
compute_becsum.o \
compute_deff.o \
compute_dip.o \
compute_rho.o \
compute_qdipol.o \
compute_qdipol_so.o \
compute_ux.o \
coset.o \
Coul_cut_2D.o \
d_matrix.o \
data_structure.o \
divide_class.o \
divide_class_so.o \
divide_et_impera.o \
rotate_wfc.o \
run_driver.o \
dynamics_module.o \
efermig.o \
efermit.o \
electrons.o \
environ_pw_module.o \
eqvect.o \
esm.o \
esm_common_mod.o \
esm_ewald_mod.o \
esm_force_mod.o \
esm_hartree_mod.o \
esm_local_mod.o \
esm_stres_mod.o \
ewald.o \
ewald_dipole.o \
extfield.o \
exx_base.o \
exx_band.o \
exx.o \
fcp_capacitance.o \
fcp_dyn_calcavg.o \
fcp_dyn_printavg.o \
fcp_dynamics.o \
fcp_hessian.o \
fcp_input.o \
fcp_module.o \
fcp_relaxation.o \
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
gcscf_input.o \
gcscf_module.o \
gen_at_dj.o \
gen_at_dy.o \
get_locals.o \
gk_sort.o \
gradcorr.o \
gweights.o \
g2_kin.o \
hs_psi.o \
hs_1psi.o \
h_epsi_her_apply.o \
h_epsi_her_set.o \
h_psi.o \
h_psi_meta.o \
hinit0.o \
hinit1.o \
hubbard.o \
init_ns.o \
init_q_aeps.o \
init_run.o \
init_us_2.o \
init_vloc.o \
input.o \
io_rho_xml.o \
irrek.o \
iweights.o \
intersite_V.o \
init_nsg.o \
nsg_adj.o \
start_k.o \
kpoint_grid.o \
lchk_tauxk.o \
ldaU.o \
make_pointlists.o \
manypw.o \
martyna_tuckerman.o \
memory_report.o \
mix_rho.o \
move_ions.o \
multable.o \
n_plane_waves.o \
new_ns.o \
new_nsb.o \
new_nsg.o \
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
plugin_check.o \
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
pw_restart_new.o \
add_qexsd_step.o \
pw_init_qexsd_input.o \
pwcom.o \
pw2blip.o \
pw2casino.o \
pw2casino_write.o \
rdiagh.o \
read_conf_from_file.o \
read_file_new.o \
realus.o \
remove_atomic_rho.o \
report_mag.o \
restart_in_electrons.o \
rho2zeta.o \
rism_module.o \
ruotaijk.o \
run_pwscf.o \
s_1psi.o \
s_psi.o \
s_psi_acc.o \
save_in_cbands.o \
save_in_electrons.o \
scale_h.o \
loc_scdm.o \
loc_scdm_k.o \
scf_mod.o \
set_kplusq.o \
set_kup_and_kdw.o \
set_occupations.o \
set_rhoc.o \
set_spin_vars.o \
set_vrs.o \
setlocal.o \
setup.o \
scissor.o \
sic.o \
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
stres_mgga.o \
stress.o \
struct_fact.o \
sum_band.o \
sumkg.o \
sumkt.o \
summary.o \
symme.o \
symm_base.o \
symmetrize_at.o \
tetra.o \
transform_becsum_so.o \
transform_becsum_nc.o \
trnvecc.o \
two_chem.o\
update_pot.o \
us_exx.o \
usnldiag.o \
v_of_rho.o \
vcsmd.o \
vcsubs.o \
vhpsi.o \
vloc_psi.o \
utils.o \
xdm_dispersion.o \
wfcinit.o \
write_ns.o \
wsweight.o \
weights.o \
ortho_wfc.o \
wannier_init.o \
wannier_check.o \
wannier_clean.o \
wannier_occ.o \
oscdft_input.o \
oscdft_context.o \
oscdft_base.o \
oscdft_enums.o \
oscdft_indices.o \
oscdft_wavefunction.o \
oscdft_wavefunction_subs.o \
oscdft_functions.o \
oscdft_occupations.o \
oscdft_wfcO.o \
oscdft_forces.o \
oscdft_forces_subs.o \
d3hess_mod.o


# GPU versions of routines
PWLIBS += \
  g_psi_mod_gpu.o \
  g_psi_gpu.o \
  h_psi_gpu.o \
  vhpsi_gpu.o \
  vloc_psi_gpu.o \
  scf_mod_gpu.o \
  usnldiag_gpu.o \
  add_vuspsi_gpu.o \
  sum_band_gpu.o \
  newd_gpu.o \
  add_paw_to_deeq_gpu.o \
  add_vhub_to_deeq_gpu.o \
  rotate_wfc_gpu.o \
  hs_1psi_gpu.o \
  hs_psi_gpu.o \
  s_1psi_gpu.o \
  utils_gpu.o \
  oscdft_functions_gpu.o



QEMODS=../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a $(BASEMODS)

TLDEPS= bindir pwlibs mods libks_solvers dftd3

all : tldeps pw.x

pw.x : $(PWOBJS) libpw.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	   $(PWOBJS) libpw.a $(QEMODS) $(QELIBS)
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

include make.depend

# DO NOT DELETE
