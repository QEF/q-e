# Makefile for 3rd derivative calculations - D3

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG).
IFLAGS=

LIBOBJS        = ../../clib/clib.a ../../iotk/src/libiotk.a

D3OBJS = \
allocate_d3.o \
allocate_pert_d3.o \
bcast_d3_input.o \
ch_psi_all2.o \
close_open.o \
d0rhod2v.o \
d2mxc.o \
d3_exc.o \
d3_init.o \
d3_readin.o \
d3_recover.o \
d3_setup.o \
d3_summary.o \
d3_symdyn.o \
d3_symdynph.o \
d3_valence.o \
d3com.o \
d3dyn_cc.o \
d3ionq.o \
d3matrix.o \
d3toten.o \
d3vrho.o \
davcio_drho2.o \
dpsi_corr.o \
dpsidpsidv.o \
dpsidvdpsi.o \
dqrhod2v.o \
drho_cc.o \
drho_drc.o \
drhod2v.o \
dvdpsi.o \
dvscf.o \
gen_dpdvp.o \
gen_dwfc.o \
incdrhoscf2.o \
openfild3.o \
print_clock_d3.o \
psymd0rho.o \
qstar_d3.o \
read_ef.o \
rotate_and_add_d3.o \
set_d3irr.o \
set_efsh.o \
set_sym_irr.o \
solve_linter_d3.o \
stop_d3.o \
sym_def1.o \
symd0rho.o \
trntnsc_3.o \
w_1gauss.o \
write_aux.o \
write_d3dyn.o \
write_igk.o \
writed3dyn_5.o

PHOBJS = ../../PHonon/PH/libph.a
LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

all : tldeps d3.x

d3.x : $(D3OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o d3.x $(D3OBJS) \
		$(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/D3/d3.x . )

tldeps :
	( cd .. ; $(MAKE) phonon || exit 1 )

clean :
	- /bin/rm -f d3.x *.o *~ *_tmp.f90 *.d *.i *.mod *.L
	- /bin/rm -f ../../bin/d3.x

include make.depend
