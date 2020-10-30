# Makefile for all_currents
# Author G. Stenuit and L. Martin-Samos

include ../../make.inc
# location of include files
#IFLAGS=-I../../include
F90FLAGS+=-fsanitize=address
LDFLAGS+=-fsanitize=address
# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../CPV/src $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../UtilXlib $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../upflib

#location of needed libraries
LIBOBJS=  ../../clib/clib.a ../../UtilXlib/libutil.a ../../KS_Solvers/DENSE/libdense.a ../../dft-d3/libdftd3qe.a ../../KS_Solvers/PPCG/libppcg.a ../../KS_Solvers/ParO/libparo.a $(QELIBS)

ALLCURRENTSOBJS = \
hartree_mod.o\
zero_mod.o\
all_currents.o \
project.o\
init_us_3.o \
init_us_1a.o \
init_reciprocal_parts_tab.o\
allocate_zero.o \
deallocate_zero.o \
ec_functionals.o \
compute_charge.o \
averages.o

TEST_CPV_TRAJ_OBJ= \
traj_object.o \
cpv_traj.o


QEMODS=   ../../Modules/libqemod.a ../../KS_Solvers/Davidson/libdavid.a ../../KS_Solvers/CG/libcg.a \
       ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a \
       ../../upflib/libupf.a
PWOBJS =  ../../LR_Modules/liblrmod.a  ../../PHonon/PH/libph.a  ../../PW/src/libpw.a ../../CPV/src/libcp.a

TLDEPS= bindir mods libs pw ph cp 

all : tldeps all_currents.x

cpv-traj-test : cpv-traj-test.x cpv_traj_test.o

cpv-traj-test.x : $(TEST_CPV_TRAJ_OBJ)
	$(LD) $(LDFLAGS) -o $@ cpv_traj_test.o $(TEST_CPV_TRAJ_OBJ) $(QEMODS)

all_currents.x : all_currents.o lib_all_currents.a $(ALLCURRENTSOBJS) $(PWOBJS) $(QEMODS)   $(TEST_CPV_TRAJ_OBJ)
	$(LD) $(LDFLAGS) -o $@ \
		all_currents.o lib_all_currents.a  $(TEST_CPV_TRAJ_OBJ) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(QEMODS) $(LIBOBJS)
	- ( cd ../../bin ; ln -fs ../energy_current/src/$@ . )

 
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_all_currents.a : $(ALLCURRENTSOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L lib_all_currents.a

maroni-sync:
	( cd ../.. ; rsync -mrvt --include="*/" --include="*Makefile" --include="*make.depend" --include="*.f90" --exclude="*" . maroni:/marconi_work/Sis19_baroni/qe-thermalcurrents2/ )

mugon-sync:
	( cd ../.. ; rsync -mrvt --include="*/" --include="*Makefile" --include="*make.depend" --include="*.f90" --exclude="*" . mugon:/scratch/rbertoss/qe-thermalcurrents2/ )

include make.depend
