# Makefile for PW

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../EE $(MOD_FLAG).

PWOBJS = \
pwscf.o \
compute_fes_grads.o \
compute_scf.o

PWLIBS = \
a2fmod.o \
add_bfield.o \
add_efield.o \
add_vuspsi.o \
addusdens.o \
addusforce.o \
addusstress.o \
allocate_bp_efield.o \
allocate_fft.o \
allocate_locpot.o \
allocate_nlpot.o \
allocate_wfc.o \
atomic_rho.o \
atomic_wfc.o \
average_pp.o \
becmod.o \
newd.o \
bp_c_phase.o \
bp_calc_btq.o \
bp_qvan3.o \
bp_strings.o \
buffers.o \
c_bands.o \
c_phase_field.o \
ccgdiagg.o \
cdiagh.o \
cdiaghg.o \
cegterg.o \
cft3.o \
cft3s.o \
checkallsym.o \
checksym.o \
clean_pw.o \
close_files.o \
commutator_Hx_psi.o \
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
davcio.o \
deriv_drhoc.o \
divide_class.o \
divide_class_so.o \
realus.o \
diropn.o \
divide.o \
divide_et_impera.o \
dqvan2.o \
drhoc.o \
dvloc_of_g.o \
dynamics_module.o \
efermig.o \
efermit.o \
electrons.o \
eqvect.o \
ewald.o \
ewald_dipole.o \
exx.o \
exx_loop.o \
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
ggen.o \
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
init_run.o \
init_us_1.o \
init_us_2.o \
init_at_1.o \
init_vloc.o \
input.o \
interpolate.o \
io_rho_xml.o \
ions.o \
irrek.o \
iweights.o \
kpoint_grid.o \
lchk_tauxk.o \
make_pointlists.o \
makov_payne.o \
martyna_tuckerman.o \
memory_report.o \
mix_pot.o \
mix_rho.o \
move_ions.o \
multable.o \
n_plane_waves.o \
new_ns.o \
new_occ.o \
ns_adj.o \
noncol.o \
non_scf.o \
offset_atom_wfc.o \
openfil.o \
orthoatwfc.o \
output_tau.o \
para.o \
paw_init.o \
paw_onecenter.o \
potinit.o \
print_clock_pw.o \
print_ks_energies.o \
punch.o \
pw_restart.o \
pwcom.o \
qvan2.o \
rcgdiagg.o \
rdiagh.o \
rdiaghg.o \
read_conf_from_file.o \
read_file.o \
read_pseudo.o \
regterg.o \
remove_atomic_rho.o \
report_mag.o \
reset_k_points.o \
restart_from_file.o \
restart_in_electrons.o \
restart_in_ions.o \
rgen.o \
rho2zeta.o \
rotate_wfc.o \
rotate_wfc_k.o \
rotate_wfc_gamma.o \
ruotaijk.o \
s_1psi.o \
s_axis_to_ca.o \
s_psi.o \
save_in_cbands.o \
save_in_electrons.o \
save_in_ions.o \
scale_h.o \
scf_mod.o \
seqopn.o \
set_fft_dim.o \
set_hubbard_l.o \
set_kplusq.o \
set_kup_and_kdw.o \
set_rhoc.o \
set_vrs.o \
setlocal.o \
setqf.o \
setup.o \
setupkpt.o \
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
stress.o \
struct_fact.o \
sum_band.o \
sumkg.o \
sumkt.o \
summary.o \
symme.o \
symm_base.o \
symmetrize_at.o \
symv.o \
tabd.o \
transform_becsum_so.o \
transform_becsum_nc.o \
trnvecc.o \
tweights.o \
update_pot.o \
usnldiag.o \
v_of_rho.o \
vcsmd.o \
vcsubs.o \
vhpsi.o \
vloc_of_g.o \
vloc_psi.o \
w0gauss.o \
w1gauss.o \
wfcinit.o \
wgauss.o \
write_ns.o \
wsweight.o \
weights.o \
ortho_wfc.o \
wannier_proj.o \
wannier_init.o \
wannier_check.o \
wannier_clean.o \
wannier_occ.o \
wannier_enrg.o 

EEOBJS=../EE/libee.a
QEMODS=../Modules/libqemod.a

TLDEPS=bindir mods libs liblapack libblas eelib

all : tldeps pw.x

pw.x : $(PWOBJS) libpw.a
	$(LD) $(LDFLAGS) -o $@ \
	   $(PWOBJS) libpw.a $(EEOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin; ln -fs ../PW/$@ . )

libpw.a : $(PWLIBS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f pw.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

include make.depend

# DO NOT DELETE
