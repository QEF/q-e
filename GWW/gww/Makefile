# Makefile for GWW
# Author: G. Stenuit and L. Martin-Samos
#
include ../../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../pw4gww

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


QEMODS = ../../KS_Solvers/libks_solvers.a $(BASEMODS)

LIBMIN= ../minpack/minpacklib.a

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS= phlibs
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../../install/tldeps

all : gww.x gww_fit.x libgww.a

gww.x : gww.o $(GWWOBJS) $(QEMODS)  $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
		gww.o $(GWWOBJS) $(QEMODS) $(LIBMIN) $(QELIBS)
	- ( cd ../../bin; ln -fs ../GWW/gww/$@ . )

gww_fit.x : gww_fit.o $(GWWOBJS) $(QEMODS) $(LIBMIN)
	$(MPIF90) $(LDFLAGS) -o $@ \
	gww_fit.o $(GWWOBJS) $(QEMODS) $(LIBMIN) $(QELIBS)
	- ( cd ../../bin; ln -fs ../GWW/gww/$@ . )

libgww.a : $(GWWOBJS)
	$(AR) $(ARFLAGS) $@ $(GWWOBJS)

clean :
	- /bin/rm -f -v gww.x gww_fit.x *.o *~ *_tmp.f90 *.d *.mod *.i *.a

include make.depend
