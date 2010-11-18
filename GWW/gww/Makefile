# Makefile for GWW
# Author: G. Stenuit and L. Martin-Samos
#
include ../../make.sys

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../pw4gww $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../flib/flib.a \
         ../../clib/clib.a ../../flib/ptools.a 

GWWOBJS = \
basic_structures.o \
compact_product.o \
create_hf.o \
create_quasi_particle.o \
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
para_gww.o \
polarization.o \
read_data_pw.o \
remainder.o \
self_energy.o \
self_energy_storage.o \
start_end.o \
times_gw.o \
vcprim.o

QEMODS = ../../Modules/libqemod.a
LIBMIN= ../minpack/minpacklib.a

TLDEPS=bindir mods libs libiotk pack

all : tldeps gww.x gww_fit.x

gww.x : gww.o $(GWWOBJS) $(LIBOBJS) $(QEMODS) $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
		gww.o $(GWWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../../bin; ln -fs ../GWW/gww/$@ . )

gww_fit.x : gww_fit.o $(GWWOBJS) $(LIBOBJS) $(QEMODS) $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
	gww_fit.o $(GWWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../../bin; ln -fs ../GWW/gww/$@ . )


tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :


clean :
	- /bin/rm -f gww.x gww_fit.x *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
