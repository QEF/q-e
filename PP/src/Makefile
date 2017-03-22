# Makefile for postprocessing (PP)

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../PW/src $(MOD_FLAG).
IFLAGS=

PPOBJS = \
addusdens1d.o \
add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
atomic_wfc_nc_proj.o \
cft.o \
chdens_module.o \
chdens_bspline.o \
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
fermisurfer_common.o \
ggen1d.o \
local_dos.o \
local_dos_mag.o \
local_dos1d.o \
openfil_pp.o\
partialdos.o \
paw_postproc.o \
plot_whole_cell.o \
projections_mod.o \
projwfc_box.o \
punch_plot.o \
smallgk.o \
stm.o \
stop_pp.o \
sym_band.o \
work_function.o \
write_p_avg.o \
write_proj.o \
write_io_header.o \
write_hamiltonians.o 

PWOBJS = ../../PW/src/libpw.a
QEOBJS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

MODULES = $(PWOBJS) $(QEOBJS)

TLDEPS= bindir mods libs pw

all : tldeps average.x bands.x dos.x epsilon.x initial_state.x fs.x \
      plan_avg.x plotband.x plotproj.x plotrho.x pmw.x pp.x projwfc.x \
      pawplot.x sumpdos.x pw2wannier90.x pw_export.x pw2gw.x \
      wannier_ham.x wannier_plot.x molecularpdos.x \
      pw2bgw.x bgw2pw.x wfck2r.x importexport_binary.x \
      fermi_velocity.x fermi_proj.x


libpp.a : $(PPOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

qexml.x: qexml.o qexml_example.o
	$(LD) $(LDFLAGS) -o $@ qexml_example.o qexml.o $(QEOBJS) $(LIBOBJS)

average.x : average.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		average.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

bands.x : bands.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		bands.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

dos.x : dos.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		dos.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pawplot.x : pawplot.o libpp.a libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pawplot.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)

plan_avg.x : plan_avg.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plan_avg.o libpp.a $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotband.x : plotband.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plotband.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotproj.x : plotproj.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plotproj.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

plotrho.x : plotrho.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		plotrho.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pmw.x : poormanwannier.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		poormanwannier.o libpp.a $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

xctest.x : xctest.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		xctest.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pp.x : postproc.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		postproc.o libpp.a $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

importexport_binary.x : importexport_binary.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		importexport_binary.o libpp.a $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

projwfc.x : projwfc.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		projwfc.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2wannier90.x : pw2wannier90.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw2wannier90.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wfck2r.x : wfck2r.o libpp.a $(MODULES) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ \
		wfck2r.o libpp.a $(MODULES) $(LIBOBJS)  $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

initial_state.x : initial_state.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		initial_state.o libpp.a $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2gw.x : pw2gw.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw2gw.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw_export.x : pw_export.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw_export.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

sumpdos.x : sumpdos.o 
	$(LD) $(LDFLAGS) -o $@ sumpdos.o $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

epsilon.x : epsilon.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ epsilon.o libpp.a $(MODULES) \
	$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wannier_ham.x : wannier_ham.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	wannier_ham.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

wannier_plot.x : wannier_plot.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
	wannier_plot.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

molecularpdos.x : molecularpdos.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		molecularpdos.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

pw2bgw.x : pw2bgw.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pw2bgw.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

bgw2pw.x : bgw2pw.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		bgw2pw.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

fs.x : fermisurface.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		fermisurface.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

fermi_velocity.x : fermi_velocity.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		fermi_velocity.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

fermi_proj.x : fermi_proj.o libpp.a $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		fermi_proj.o libpp.a $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PP/src/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *_tmp.f90 *.d *.mod *.i *.L libpp.a

include make.depend
