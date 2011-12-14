# Makefile for ballistic conductance (PWCOND)

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG).

PPOBJS = \
../../PP/xsf.o

PWCONDOBJS = \
allocate_cond.o \
bessj.o \
compbs.o \
compbs_2.o \
condcom.o \
condmain.o \
cond_out.o \
cond_restart.o\
do_cond.o \
eigenchnl.o \
form_zk.o \
four.o \
free_mem.o \
gep_x.o \
gramsh.o \
hev_ab.o \
init_cond.o \
init_orbitals.o \
init_gper.o \
integrals.o \
jbloch.o \
kbloch.o \
local.o \
local_set.o \
openfil_cond.o \
poten.o \
print_clock_pwcond.o \
realus_scatt.o \
rotproc.o \
save_cond.o \
scat_states_plot.o \
scatter_forw.o \
summary_band.o \
summary_tran.o \
sunitary.o \
transmit.o 

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a

LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

TLDEPS= bindir mods libs pw pp

all : tldeps pwcond.x

pwcond.x : $(PWCONDOBJS) $(PWOBJS) $(PPOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o pwcond.x $(PWCONDOBJS) $(PPOBJS) \
		$(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PWCOND/src/pwcond.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :
 
clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.a *.L
	- /bin/rm -f ../../bin/pwcond.x

include make.depend
