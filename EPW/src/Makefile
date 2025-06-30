# Makefile for EPW

include ../../make.inc

#
# use recursive definitions in Makefile to modify the PATH for modules
#

MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../external/wannier90/src/obj \
          $(MOD_FLAG)../../dft-d3

#
# The above is in alphabetical order except when order matter during compilation
#

EPWOBJS = \
input.o \
ep_constants.o \
global_var.o \
io/io_var.o \
bloch2wannier.o \
utilities/low_lvl.o \
utilities/parallelism.o \
longrange.o \
io/io.o \
io/io_transport.o \
io/io_selfen.o \
wigner.o \
wannier2bloch.o \
wannier2bloch_opt.o \
cumulant.o \
close.o \
printing.o \
utilities/kfold.o \
utilities/bzgrid.o \
io/io_supercond.o \
supercond_vertex.o \
utilities/utilities.o \
supercond_common.o \
supercond.o \
supercond_iso.o \
supercond_aniso.o \
supercond_coul.o \
utilities/sparse_ir.o \
io/io_sparse_ir.o \
stop.o \
transport_mag.o \
utilities/symmetry.o \
transport_legacy.o \
transport.o \
pw2wan.o \
wannierization.o \
indabs.o \
io/io_indabs.o\
utilities/bcast_input.o \
io/io_ahc.o \
dvqpsi.o \
ep_coarse.o \
ep_coarse_unfolding.o \
supercond_driver.o \
selfen.o \
spectral.o \
wannier.o \
use_wannier.o \
init.o \
readin.o \
setups.o \
summaries.o \
wann_common.o \
polaron.o \
qdabs.o\
wfpt.o \
exphonon.o \
expolaron.o \
utilities/screening.o \
tdbe_driver.o \
tdbe_mod.o \
tdbe_common.o \
phph.o \
utilities/screening.o \
#default : epw

PHOBJS = ../../PHonon/PH/libph.a
W90LIB = ../../external/wannier90/libwannier.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

TLDEPS= pw ph pp

all : pw ph pp wannier libepw.a epw.x nscf2supercond.x

libepw.a : $(EPWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

epw.x : epw.o libepw.a $(PHOBJS) $(LRMODS) $(PWOBJS) $(W90LIB) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	epw.o  libepw.a	$(PHOBJS) $(LRMODS) $(W90LIB) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../bin ; ln -fs ../src/epw.x . )

nscf2supercond.x : utilities/nscf2supercond.o $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	utilities/nscf2supercond.o ../../PP/src/libpp.a $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../bin ; ln -fs ../src/$@ . )
	- ( cd ../../bin ; ln -fs ../EPW/src/$@ . )

pw :
	cd ../../ ; make pw

ph :
	cd ../../ ; make ph

pp :
	cd ../../ ; make pp

wannier :
	cd ../../ ; make w90

ld1 :
	cd ../../ ; make ld1

clean :
	- /bin/rm -f  *.o *~ *.d *.mod *.i io/*.mod io/*.o utilities/*.mod utilities/*.o libepw.a liblr.a

include make.depend
