# Makefile for pw4gww (PW4GWW)
# Author G. Stenuit and L. Martin-Samos

include ../../make.sys

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../flib/flib.a \
         ../../clib/clib.a ../../flib/ptools.a

PW4GWWOBJS = \
allocate_wannier.o \
calculate_wing.o \
cprim.o \
dft_exchange.o \
dirdel.o \
distance_gw.o \
energies_u.o \
energies_xc.o \
go_wannier_product.o \
gram_schmidt_pwannier.o \
local_wannier.o \
matrix_wannier_gamma.o \
memory_limit.o \
produce_wannier_gamma.o \
product_wannier.o \
product_wannier_para.o \
product_wannier_para_c.o \
openfil_pw4gww.o \
stop_pp.o \
real_wfc.o \
rotate_wannier.o \
set_wannier_P.o \
start_pw4gww.o \
ultra_external.o \
ultralocalization.o \
ultralocalization_para.o \
wannier.o \
wannier_gw.o \
wannier_pmat_terms_ggrid.o \
wannier_uterms.o \
wannier_uterms_c.o \
wannier_valence_terms.o \
wfc_real.o \
write_vpot_matrix.o \
write_wannier_matrix.o \
write_wannier_matrix_c.o \
write_wannier_products.o 

QEMODS = ../../Modules/libqemod.a
PWOBJS = ../../PW/src/libpw.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS= bindir libs pw

all : tldeps pw4gww.x

pw4gww.x : pw4gww.o libpw4gww.a $(PW4GWWOBJS) $(PWOBJS) $(QEMODS) $(LIBMIN) 
	$(LD) $(LDFLAGS) -o $@ \
		pw4gww.o libpw4gww.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../GWW/pw4gww/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libpw4gww.a : $(PW4GWWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libpw4gww.a

include make.depend
