# Makefile for EPW

include ../../make.inc
include make.libs

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
##### io_scattering.o \

EPWOBJS = \
epwcom.o \
constants_epw.o \
elph2.o \
io_var.o \
low_lvl.o \
division.o \
rigid_epw.o \
io_epw.o \
io_transport.o \
io_selfen.o \
wigner.o \
wan2bloch.o \
cum_mod.o \
close_epw.o \
poolgathering.o \
printing.o \
kfold.o \
dynmat_asr.o \
io_eliashberg.o \
utilities.o \
eliashbergcom.o \
supercond.o \
supercond_iso.o \
supercond_aniso.o \
stop_epw.o \
grid.o \
bfield.o \
rotate.o \
rotate_wavefunction.o \
transport.o \
transport_iter.o \
pw2wan2epw.o \
wannierization.o \
indabs.o \
io_indabs.o\
bcast_epw_input.o \
bloch2wan.o \
read_ahc_files.o \
dvqpsi.o \
eliashberg.o \
elphel2_shuffle.o \
elphon_shuffle.o \
elphon_shuffle_wrap.o \
selfen.o \
selfen_ahc.o \
spectral_func.o \
ephwann_shuffle.o \
epw_init.o \
epw_readin.o \
epw_setup.o \
epw_summary.o \
loadumat.o \
wannierEPW.o \
loadbm.o \
polaron.o \
qdabs.o\
wfpt_mod.o \
#default : epw


PHOBJS = ../../PHonon/PH/libph.a
W90LIB = ../../external/wannier90/libwannier.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)
#SCALAPACK_LIBS = -lmkl_blacs_intelmpi_lp64 -lmkl_scalapack_lp64

TLDEPS= pw ph pp

all :   pw ph pp wannier libepw.a epw.x

libepw.a : $(EPWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

epw.x : epw.o libepw.a $(PHOBJS) $(LRMODS) $(PWOBJS) $(W90LIB) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	epw.o  libepw.a	$(PHOBJS) $(LRMODS) $(W90LIB) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../bin ; ln -fs ../src/epw.x . )

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
	- /bin/rm -f  *.o *~ *.d *.mod *.i libepw.a liblr.a

include make.depend
