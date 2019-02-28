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
          $(MOD_FLAG)../../wannier90-3.0.0/src/obj \
          $(MOD_FLAG)../../dft-d3

#
# The above is in alphabetical order except when order matter during compilation
#

EPWOBJS = \
epwcom.o \
constants_epw.o \
io_epw.o \
elph2.o \
division.o \
eliashbergcom.o \
superconductivity.o \
superconductivity_aniso.o \
superconductivity_iso.o \
io_eliashberg.o \
io_scattering.o \
transportcom.o \
close_epw.o \
printing.o \
transport.o \
transport_iter.o \
wigner.o \
a2f.o \
adddvscf2.o \
allocate_epwq.o \
bcast_epw_input.o \
bloch2wan.o \
broyden.o \
createkmap.o \
dvanqq2.o \
dvqpsi_us3.o \
dvqpsi_us_only3.o \
eliashberg.o \
elphel2_shuffle.o \
elphon_shuffle.o \
elphon_shuffle_wrap.o \
ephwann_shuffle.o \
ephwann_shuffle_mem.o \
epw_init.o \
epw_readin.o \
epw_setup.o \
epw_summary.o \
fermiwindow.o \
gmap_sym.o \
indabs.o \
io_dyn_mat2.o \
kfold.o \
kpoint_grid_epw.o \
ktokpmq.o \
loadkmesh.o \
loadqmesh.o \
loadumat.o \
load_rebal.o \
nesting_fn.o \
newdq2.o \
openfilepw.o \
rgd_blk_epw_fine_mem.o \
plot_band.o \
poolgather.o \
print_ibte.o \
readdvscf.o \
readgmap.o \
readmat_shuffle2.o \
readwfc.o \
rigid_epw.o \
rotate_eigenm.o \
rotate_epmat.o \
selfen_elec.o \
selfen_phon.o \
selfen_pl.o \
set_ndnmbr.o \
setphases_wrap.o \
sort.o \
spectral_cumulant.o \
spectral_func.o \
spectral_func_ph.o \
spectral_func_pl.o \
star_q2.o \
stop_epw.o \
wannierEPW.o \
wannierize.o \
pw2wan90epw.o \
wan2bloch.o \
system_mem_usage.o

# setphases.o \ : SP: We keep it in case its usefull in the future but now depreciated.

#default : epw


PHOBJS = ../../PHonon/PH/libph.a
PWOBJS = ../../PW/src/libpw.a
W90LIB = ../../wannier90-3.0.0/libwannier.a 
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a 
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
         ../../FFTXlib/libqefft.a ../../dft-d3/libdftd3qe.a
LIBOBJS =../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../clib/clib.a ../../iotk/src/libiotk.a 

TLDEPS= pw ph

all :   pw ph wannier pp upf libepw.a epw.x 

libepw.a : $(EPWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

epw.x : epw.o libepw.a $(PHOBJS) $(LRMODS) $(PWOBJS) $(W90LIB) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	epw.o  libepw.a	$(PHOBJS) $(LRMODS) $(W90LIB) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) 
	- ( cd ../bin ; ln -fs ../src/epw.x . )

pw : 
	cd ../../ ; make pw

ph : 
	cd ../../ ; make ph

wannier :
	cd ../../ ; make w90 ; cd wannier90-3.0.0/ ; make lib

pp : 
	cd ../../ ; make pp

ld1 :
	cd ../../ ; make ld1

upf :
	cd ../../ ; make upf

clean :
	- /bin/rm -f  *.o *~ *.d *.mod *.i libepw.a liblr.a

include make.depend
