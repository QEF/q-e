# Makefile for KCW. 

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG).


KI_OBJS= \
wann2kcw.o\
kcw_screen.o\
kcw_ham.o\
kcw_comm.o \
print_clock_kcw.o \
kcw_readin.o \
kcw_setup.o \
kcw_setup_screen.o \
kcw_setup_ham.o \
bcast_kcw_input.o \
input_summary.o \
read_wannier.o \
rotate_ks.o \
apply_u_matrix.o \
ks_hamiltonian.o \
screen_coeff.o \
compute_map_ikq_single.o \
calculate_phase.o \
structure_factor.o \
rho_of_q.o \
kcw_prepare_q.o \
bare_pot.o \
solve_linter_koop_mod.o \
kcw_initialize_ph.o \
kcw_openfilq.o \
kcw_q_setup.o \
clean_pw_kcw.o \
close_kcw.o \
kcw_R_points.o\
ham_R0_2nd.o\
read_alpha.o\
kcw_io_new.o\
bcast_wfc.o\
interpolation.o\
convert_kpts_names.o\
write_hr_to_file.o\
alpha_corr.o \
koopmans_ham.o \
full_ham.o \
self_hartree.o\
group_orbitals.o\
kcw_run_nscf.o\
kcw_init_q.o\
kcw_allocate_q.o\
kcw_deallocate_q.o\
setup_coulomb_exx.o\
coulomb.o
#restore_symm.o 
#check_density.o 
#setup_coulomb.o\

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../upflib/libupf.a ../../XClib/xc_lib.a \
         ../../FFTXlib/src/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a \
          ../..//MBD/libmbd.a

TLDEPS= pwlibs lrmods pp w90 #Dependency. To be compiled first

all : tldeps libs-kcw kcw.x clean_partial

libs-kcw : libkcw.a

libkcw.a : $(KI_OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

kcw.x:  kcw.o libkcw.a $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ kcw.o libkcw.a\
        $(LRMODS) $(PWOBJS) $(QEMODS)  \
        $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../KCW/src/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *__genmod.f90 *__genmod.mod ../../bin/kcw.x

clean_partial :
	- /bin/rm -f *__genmod.f90 *__genmod.mod

include make.depend
# DO NOT DELETE
