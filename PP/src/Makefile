# Makefile for postprocessing (PP)

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG).

PPOBJS = \
addusdens1d.o \
add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
atomic_wfc_nc_proj.o \
cft.o \
cgracsc.o \
compute_ppsi.o \
compute_sigma_avg.o \
cube.o \
spinsym.o \
dosg.o \
do_initial_state.o \
do_shift_ew.o \
d_matrix_nc.o \
d_matrix_so.o \
elf.o \
ggen1d.o \
local_dos.o \
local_dos_mag.o \
local_dos1d.o \
openfil_pp.o\
paw_postproc.o \
plot_whole_cell.o \
punch_plot.o \
smallgk.o \
stm.o \
stop_pp.o \
sym_band.o \
work_function.o \
write_p_avg.o \
xsf.o \
wannier_hamiltonians.o 

PWOBJS = ../../PW/src/libpw.a
QEOBJS = ../../Modules/libqemod.a

MODULES = $(PWOBJS) $(QEOBJS)

TLDEPS= bindir mods libs pw

LIBOBJS = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

all : tldeps average.x bands.x dos.x epsilon.x initial_state.x \
      plan_avg.x plotband.x plotproj.x plotrho.x pmw.x pp.x projwfc.x \
      pawplot.x sumpdos.x pw2wannier90.x pw_export.x pw2gw.x \
      wannier_ham.x wannier_plot.x \
      pw2bgw.x bgw2pw.x

qexml.x: qexml.o qexml_example.o
	$(LD) $(LDFLAGS) -o $@ qexml_example.o qexml.o $(QEOBJS) $(LIBOBJS)

average.x : average.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		average.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

bands.x : bands.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		bands.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

dos.x : dos.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		dos.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pawplot.x : pawplot.o $(PPOBJS) $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pawplot.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)

plan_avg.x : plan_avg.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plan_avg.o $(PPOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotband.x : plotband.o $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plotband.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotproj.x : plotproj.o $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plotproj.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotrho.x : plotrho.o $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plotrho.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pmw.x : poormanwannier.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		poormanwannier.o $(PPOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

xctest.x : xctest.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		xctest.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pp.x : postproc.o chdens.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		postproc.o chdens.o $(PPOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

projwfc.x : projwfc.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		projwfc.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2wannier90.x : pw2wannier90.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw2wannier90.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

initial_state.x : initial_state.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		initial_state.o $(PPOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2gw.x : pw2gw.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw2gw.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw_export.x : pw_export.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw_export.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

sumpdos.x : sumpdos.o 
	$(LD) $(LDFLAGS) -o $@ sumpdos.o
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

epsilon.x : epsilon.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ epsilon.o $(PPOBJS) $(MODULES) \
	$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wannier_ham.x : wannier_ham.o wannier_umatrix_subroutines.o wannier_u_matrix.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	wannier_ham.o wannier_umatrix_subroutines.o wannier_u_matrix.o \
	$(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )
        
wannier_plot.x : wannier_plot.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	wannier_plot.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2bgw.x : pw2bgw.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw2bgw.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

bgw2pw.x : bgw2pw.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		bgw2pw.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L

include make.depend
