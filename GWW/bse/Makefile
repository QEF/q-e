# Makefile for BSE
#
 include ../../make.inc

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PHonon/Gamma/ $(MOD_FLAG)../../PHonon/PH/ \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../pw4gww \
          $(MOD_FLAG)../gww $(MOD_FLAG).

BSEOBJS = \
openfil_bse.o \
start_bse.o \
bse_wannier.o\
transitions.o\
read_export.o\
bse_basic_structure.o\
exciton.o\
diago_exc.o\
exchange_exc.o\
rotate_wannier.o\
write_wannier_matrix.o\
direct_v_exc.o\
direct_w_exc.o\
exc_h_a.o\
sdescent.o\
check_basis.o\
conjgrad.o\
conj_grad_stuff.o\
find_eig.o\
plot_excwfn.o\
tspace.o\
absorption.o\
spectrum.o\
absorption.o\
cgsolve.o\
h_h.o\
dvpsi_bse.o\
lanczos.o\
print_spectrum.o\
print_bse.o\
contract_w_exc.o\
qpe_exc.o\
qpcorrections.o
#zvscal_test.o\


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a
PW4GWWOBJ = ../pw4gww/fft_custom.o ../pw4gww/stop_pp.o ../pw4gww/mp_wave_parallel.o
GWWOBJ = ../gww/libgww.a ../minpack/minpacklib.a
#LIBPWPH =  ../../PHonon/Gamma/libcg.a
LIBPWPH =
#PPOBJ= ../../PP/src/libpp.a
PPOBJ=

TLDEPS= bindir libs pw ph mods 

all : tldeps bse_main.x

bse_main.x : bse_main.o libbse.a  $(BSEOBJS) $(PWOBJS) $(QEMODS) $(LIBPWPH)  $(GWWOBJ) $(PPOBJ)
	$(LD) $(LDFLAGS) -o $@ \
		bse_main.o libbse.a $(PWOBJS) $(EEOBJS) $(LIBPWPH) $(QEMODS) $(LIBOBJS) $(LIBS)  $(PW4GWWOBJ) $(GWWOBJ) $(PPOBJ)
	- ( cd ../../bin ; ln -fs ../GWW/bse/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libbse.a : $(BSEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libbse.a

include make.depend



