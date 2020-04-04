# Makefile for all_currents
# Author G. Stenuit and L. Martin-Samos

include ../make.inc
# location of include files
#IFLAGS=-I../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG)../LAXlib \
          $(MOD_FLAG)../PW/src $(MOD_FLAG)../CPV/src $(MOD_FLAG)../PHonon/PH \
          $(MOD_FLAG)../UtilXlib $(MOD_FLAG)../FFTXlib $(MOD_FLAG)../LR_Modules 

#location of needed libraries
LIBOBJS= ../iotk/src/libiotk.a ../clib/clib.a ../UtilXlib/libutil.a 

ALLCURRENTSOBJS = \
hartree_mod.o\
zero_mod.o\
diropn_due.o\
all_currents.o \
openfil_all_currents.o\
stop_pp.o\
routine_hartree.o\
project.o\
orthogonalize2.o \
init_us_3.o \
init_us_1a.o \
add_nc_curr.o \
stop_pp.o\
routine_zero.o\
init_reciprocal_parts_tab.o\
allocate_zero.o \
deallocate_zero.o \
aris_functionals.o 



QEMODS=   ../Modules/libqemod.a ../KS_Solvers/Davidson/libdavid.a ../KS_Solvers/CG/libcg.a \
       ../FFTXlib/libqefft.a ../LAXlib/libqela.a ../UtilXlib/libutil.a 
PWOBJS =  ../LR_Modules/liblrmod.a  ../PHonon/PH/libph.a  ../PW/src/libpw.a ../CPV/src/libcp.a

TLDEPS= bindir mods libs pw ph cp 

all : tldeps all_currents.x

all_currents.x : all_currents.o lib_all_currents.a $(ALLCURRENTSOBJS) $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ \
		all_currents.o lib_all_currents.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../all_currents/$@ . )

# ( cd ../../bin ; ln -fs ../GWW/pw4gww/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

lib_all_currents.a : $(ALLCURRENTSOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L lib_all_currents.a

maroni-sync:
	rsync Makefile make.depend *.f90 maroni:/marconi_work/Sis19_baroni/qe-thermalcurrents2/all_currents

include make.depend
