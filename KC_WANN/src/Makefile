# Makefile for KC_WANN. 

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG).


KI_OBJS= \
kc_wann_com.o \
print_clock_kc.o \
kc_readin.o \
kc_setup.o \
kc_setup_screen.o \
kc_setup_ham.o \
bcast_kc_input.o \
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
kc_prepare_q.o \
bare_pot.o \
solve_linter_koop_mod.o \
kc_initialize_ph.o \
kc_openfilq.o \
kc_phq_setup.o \
clean_pw_kc.o \
close_kc.o \
kc_R_points.o\
ham_R0_2nd.o\
hamilt.o\
read_alpha.o\
setup_coulomb.o\
setup_coulomb_exx.o\
kc_io_new.o\
bcast_wfc.o\
interpolation.o\
convert_kpts_names.o\
write_hr_to_file.o\
alpha_corr.o \
koopmans_ham.o \
full_ham.o \
self_hartree.o\
group_orbitals.o
#restore_symm.o \
#check_density.o 



LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
PHOBJS = ../../PHonon/PH/libph.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../upflib/libupf.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a


TLDEPS= pw lrmods phlibs pp w90 #Dependency. To be compiled first

#all : qehacks tldeps odd_alpha.x ki.x 
all : tldeps wann2kc.x kc_screen.x kc_ham.x clean_partial

wann2kc.x:  wann2kc.o $(KI_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ wann2kc.o \
        $(KI_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)  \
        $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../KC_WANN/src/$@ . )

kc_screen.x :  kc_screen.o $(KI_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ kc_screen.o \
	$(KI_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)  \
	$(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../KC_WANN/src/$@ . )

kc_ham.x :  kc_ham.o $(COMMON_OBJS) $(KI_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ kc_ham.o \
	$(KI_OBJS) $(PHOBJS) $(LRMODS) $(PHOBJS) $(PWOBJS) $(QEMODS)  \
	$(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../KC_WANN/src/$@ . )


tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *__genmod.f90 *__genmod.mod

clean_partial :
	- /bin/rm -f *__genmod.f90 *__genmod.mod

include make.depend
# DO NOT DELETE
