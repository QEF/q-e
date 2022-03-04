#
# Copyright (C) 2001-2021 Quantum ESPRESSO group
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
#
#-----------------------------------------------------------------------
# Makefile for QEHeat (all_currents.x)

include ../../make.inc
# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../LR_Modules



TEST_CPV_TRAJ_OBJ= \
traj_object.o \
cpv_traj.o

ALLCURRENTSOBJS = \
test_h_psi_s_psi_commutator_Hx_psi.o\
hartree_xc_mod.o\
kohn_sham_mod.o\
zero_mod.o\
ionic_mod.o\
project.o\
init_us_3.o \
ec_functionals.o \
compute_charge.o \
averages.o\
scf_result.o\
$(TEST_CPV_TRAJ_OBJ)


QEMODS=../../PHonon/PH/libph.a ../../LR_Modules/liblrmod.a \
       ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a \
      $(BASEMODS)

TLDEPS= pw ph

all : tldeps all_currents.x

cpv-traj-test : cpv-traj-test.x cpv_traj_test.o

cpv-traj-test.x : $(TEST_CPV_TRAJ_OBJ)
	$(LD) $(LDFLAGS) -o $@ cpv_traj_test.o $(TEST_CPV_TRAJ_OBJ) $(QEMODS)

all_currents.x : all_currents.o lib_all_currents.a $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
		all_currents.o lib_all_currents.a $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../QEHeat/src/$@ . )


tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_all_currents.a : $(ALLCURRENTSOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L lib_all_currents.a

include make.depend
