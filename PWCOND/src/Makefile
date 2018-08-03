# Makefile for ballistic conductance (PWCOND)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG) ../../KS_Solvers/Davidson \
          $(MOD_FLAG)../../KS_Solvers/CG \
          $(MOD_FLAG)../../KS_Solvers/PPCG \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3

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
plus_u_setup.o \
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
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/Davidson/libdavid.a \
         ../../KS_Solvers/CG/libcg.a ../../KS_Solvers/PPCG/libppcg.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a \
	 ../../dft-d3/libdftd3qe.a

TLDEPS= bindir mods libs pw pp

all : tldeps pwcond.x

pwcond.x : $(PWCONDOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o pwcond.x $(PWCONDOBJS) \
		$(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PWCOND/src/pwcond.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.a *.L
	- /bin/rm -f ../../bin/pwcond.x

include make.depend
