# Makefile for GWW
# Author: G. Stenuit and L. Martin-Samos
#
include ../../make.inc

# location of include files
IFLAGS=

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../pw4gww $(MOD_FLAG).

GWWOBJS = \
basic_structures.o \
compact_product.o \
contour.o \
create_hf.o \
create_quasi_particle.o \
create_quasi_particle_off.o \
do_contour.o \
do_polarization_lanczos.o \
 do_self_lanczos.o \
do_self_lanczos_time.o \
do_self_lanczos_full.o \
energies_gww.o \
expansion.o \
fft_gw.o \
fit_multipole.o \
fit_polynomial.o \
go_dressed_w.o \
go_exchange.o \
go_fft.o \
go_green.o \
go_polarization.o \
green_function.o \
gv_time.o \
input_gw.o \
lanczos_polarization.o \
para_gww.o \
polarization.o \
read_data_pw.o \
self_energy.o \
self_energy_storage.o \
start_end.o \
times_gw.o \
vcprim.o 


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS=bindir mods libs libiotk 

all : tldeps gww.x gww_fit.x libgww.a

gww.x : gww.o $(GWWOBJS) $(LIBOBJS) $(QEMODS)  $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
		gww.o $(GWWOBJS) $(QEMODS) $(LIBOBJS) $(LIBMIN) $(LIBS)
	- ( cd ../../bin; ln -fs ../GWW/gww/$@ . )

gww_fit.x : gww_fit.o $(GWWOBJS) $(LIBOBJS) $(QEMODS) $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
	gww_fit.o $(GWWOBJS) $(QEMODS) $(LIBOBJS) $(LIBMIN) $(LIBS)
	- ( cd ../../bin; ln -fs ../GWW/gww/$@ . )

libgww.a : $(GWWOBJS)
	$(AR) $(ARFLAGS) $@ $(GWWOBJS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f -v gww.x gww_fit.x *.o *~ *_tmp.f90 *.d *.mod *.i *.a

include make.depend

# DO NOT DELETE
