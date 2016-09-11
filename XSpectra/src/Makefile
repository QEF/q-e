# Makefile for X-ray spectra calculation (XSPECTRA)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG).
IFLAGS=

XOBJS = \
./gaunt_mod.o \
./xspectra_mod.o \
./radin_mod.o \
./mygetK.o \
./ipoolscatter.o \
./lr_sm1_psi.o \
./orthoUatwfc_k.o \
./read_k_points.o \
./reset_k_points_and_reinit.o \
./select_nl_init.o \
./lanczos.o \
./plot_xanes_cross_sections.o \
./xanes_dipole_general_edge.o\
./banner_xspectra.o \
./set_xspectra_namelists_defaults.o \
./read_input_and_bcast.o \
./stdout_routines.o \
./xanes_dipole.o \
./check_orthogonality_k_epsilon.o \
./assign_paw_radii_to_species.o \
./io_routines.o \
./xanes_quadrupole.o 



GIPAWOBJS=./paw_gipaw.o \
./gipaw_module.o \
./init_gipaw_1.o \
./init_gipaw_2.o

MANIP_XS_OBJ=./gaunt_mod.o 
             
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a 

TLDEPS = bindir libs pw

all : tldeps xspectra.x spectra_correction.x molecularnexafs.x

xspectra.x :   xspectra.o $(XOBJS) $(GIPAWOBJS) $(QEMODS) $(PWOBJS)
	$(MPIF90) -o $@  xspectra.o $(XOBJS) $(GIPAWOBJS) \
		$(PWOBJS) $(QEMODS) $(LDFLAGS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../XSpectra/src/$@ . )

spectra_correction.x : spectra_correction.o $(MANIP_XS_OBJ) $(QEMODS) $(PWOBJS)
	$(MPIF90) -o $@  spectra_correction.o $(MANIP_XS_OBJ) $(PWOBJS) $(QEMODS) \
                $(LDFLAGS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../XSpectra/src/$@ . )

molecularnexafs.x : molecularnexafs.o 
	$(MPIF90) -o $@  molecularnexafs.o $(LDFLAGS)
	- ( cd ../../bin ; ln -fs ../XSpectra/src/$@ . )

tldeps : 
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/xspectra.x

include make.depend
# DO NOT DELETE
