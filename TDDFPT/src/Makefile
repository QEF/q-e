# Makefile for TDDFPT

include ../../make.inc

IFLAGS         = -I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib \
          $(MOD_FLAG) . $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../LR_Modules 

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a
LRMODS = ../../LR_Modules/liblrmod.a

LROBJS= \
lr_variables.o \
lr_charg_resp.o \
bcast_lr_input.o \
lr_readin.o \
lr_alloc_init.o \
lr_calc_dens.o \
lr_dot.o \
lr_dealloc.o \
lr_ortho.o \
lr_read_wf.o \
lr_normalise.o \
lr_lanczos.o \
lr_apply_liouvillian.o \
lr_dv_setup.o \
lr_solve_e.o \
lr_dvpsi_e.o \
stop_lr.o \
lr_read_d0psi.o \
lr_restart.o \
lr_write_restart.o \
print_clock_lr.o \
sd0psi.o \
lr_set_boxes_density.o \
lr_init_nfo.o \
lr_exx_kernel.o \
lr_calc_dens_eels.o \
lr_calc_dens_eels_nc.o \
lr_apply_liouvillian_eels.o \
lr_addusddens.o \
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
lr_dav_debug.o

TLDEPS=bindir libs pw-lib libfft libla lrmods

all : tldeps libtddfpt.a turbo_lanczos.x turbo_davidson.x turbo_eels.x

libtddfpt.a : $(LROBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

turbo_lanczos.x : lr_main.o libtddfpt.a $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		lr_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_davidson.x : lr_dav_main.o libtddfpt.a $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		lr_dav_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_eels.x : lr_eels_main.o libtddfpt.a $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
                lr_eels_main.o libtddfpt.a $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f ../bin/turbo_lanczos.x ../bin/turbo_davidson.x ../bin/turbo_eels.x *.x *.o *~ *.d *.mod *.i work.pc *.a *_tmp.f90

include make.depend

# DO NOT DELETE
