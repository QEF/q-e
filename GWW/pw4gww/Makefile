# Makefile for pw4gww (PW4GWW)
# Author G. Stenuit and L. Martin-Samos

include ../../make.inc

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG).

PW4GWWOBJS = \
allocate_wannier.o \
produce_wannier_gamma.o \
start_pw4gww.o \
stop_pp.o \
openfil_pw4gww.o \
dft_exchange.o  \
wfc_real.o \
full.o \
energies_xc.o  \
wannier.o  \
write_wannier_matrix.o  \
rotate_wannier.o  \
self_lanczos.o \
pola_lanczos.o \
matrix_wannier_gamma.o \
calculate_wing.o  \
fake_conduction.o \
o_1psi.o  \
o_bands.o  \
optimal.o \
pola_partial.o \
semicore.o \
semicore_read.o \
v_basis.o  \
wannier_uterms.o \
write_vpot_matrix.o  \
diago_cg.o  \
o_rinitcgg.o \
diago_cg_g.o \
contour_terms.o \
fft_custom.o  \
exchange_custom.o \
mp_wave_parallel.o \
wannier_bse.o\
hpsi_pw4gww.o \
cgsolve_all_gamma.o \
realus.o 

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS= bindir libs pw

all : tldeps pw4gww.x

pw4gww.x : pw4gww.o libpw4gww.a  $(PW4GWWOBJS) $(PWOBJS) $(QEMODS)  $(LIBOBJS) $(LIBMIN) 
	$(LD) $(LDFLAGS) -o $@ \
	pw4gww.o libpw4gww.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBMIN) $(LIBS)
	- ( cd ../../bin ; ln -fs ../GWW/pw4gww/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

libpw4gww.a : $(PW4GWWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -fv *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libpw4gww.a

include make.depend
# DO NOT DELETE
