# Makefile for TDDFPT

include ../../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules

QEMODS = $(BASEMODS)
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
LRMODS = ../../LR_Modules/liblrmod.a

LROBJS= \
environ_td_module.o \
lr_variables.o \
lr_charg_resp.o \
bcast_lr_input.o \
lr_readin.o \
lr_alloc_init.o \
lr_calc_dens.o \
lr_dealloc.o \
lr_ortho.o \
lr_read_wf.o \
lr_normalise.o \
lr_lanczos.o \
lr_apply_liouvillian.o \
lr_dv_setup.o \
lr_compute_intq.o \
lr_set_intq_nc.o \
lr_mix_pot.o \
paw_add_symmetry.o \
lr_transform_intq_nc.o \
lr_transform_intq_so.o  \
dveqpsi_us.o \
dveqpsi_us_only.o \
lr_solve_e.o \
lr_dvpsi_e.o \
stop_lr.o \
lr_read_d0psi.o \
lr_restart.o \
lr_write_restart.o \
print_clock_lr.o \
sd0psi.o \
lr_set_boxes_density.o \
lr_magnons_routines.o\
lr_init_nfo.o \
lr_exx_kernel.o \
lr_calc_dens_eels.o \
lr_calc_dens_eels_nc.o \
lr_apply_liouvillian_eels.o \
lr_addus_dvpsi.o \
lr_run_nscf.o \
lr_setup_nscf.o \
lr_sym_eels.o \
lr_smallgq.o \
lr_psym_eels.o \
lr_summary.o \
lr_dvpsi_eels.o \
lr_dav_variables.o\
lr_dav_routines.o\
lr_us.o\
lr_test_restart.o\
lr_dav_debug.o\
plugin_tddfpt_potential.o\
linear_solvers.o\
orthogonalize_omega.o\
lr_sternheimer.o\
lr_apply_liouvillian_magnons.o\
lr_dvpsi_magnons.o\
lr_Opsi_magnons.o\
lr_calc_dens_magnons.o

TLDEPS=lrmods

all : tldeps libtddfpt.a turbo_lanczos.x turbo_davidson.x turbo_eels.x turbo_magnon.x turbo_spectrum.x

libtddfpt.a : $(LROBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

turbo_lanczos.x : lr_main.o libtddfpt.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
		lr_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_davidson.x : lr_dav_main.o libtddfpt.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
		lr_dav_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_eels.x : lr_eels_main.o libtddfpt.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
                lr_eels_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_magnon.x : lr_magnons_main.o libtddfpt.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	        lr_magnons_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_spectrum.x : turbo_spectrum.o libtddfpt.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
		turbo_spectrum.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f ../bin/turbo_lanczos.x ../bin/turbo_davidson.x ../bin/turbo_eels.x \
		../bin/turbo_magnons.x ../bin/turbo_spectrum.x *.x *.o *~ *.d *.mod *.i work.pc *.a *_tmp.f90

include make.depend

# DO NOT DELETE
