# Makefile for PW
#Modified for v4.0 by obm

include ../../make.sys
IFLAGS         = -I../../include
LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a 


# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG) . $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../PHonon/PH

QEMODS = ../../Modules/libqemod.a
PWOBJS = ../../PW/src/libpw.a
PHOBJS=../../PHonon/PH/libph.a

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
lr_sm1_psi.o \
stop_lr.o \
lr_read_d0psi.o \
lr_restart.o \
lr_write_restart.o \
print_clock_lr.o \
sd0psi.o \
lr_set_boxes_density.o \
lr_init_nfo.o \
lr_exx_kernel.o \
lr_dav_variables.o\
lr_dav_routines.o\
lr_us.o\
lr_dav_debug.o

TLDEPS=bindir libs pw ph

all :  tldeps libtddfpt.a turbo_lanczos.x turbo_davidson.x

libtddfpt.a : $(LROBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

turbo_lanczos.x : lr_main.o libtddfpt.a $(PHOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		lr_main.o libtddfpt.a $(PHOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

turbo_davidson.x : lr_dav_main.o libtddfpt.a $(PHOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		lr_dav_main.o libtddfpt.a $(PHOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/src/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../src/$@ . ); fi

dene.x : obm_tests.o 
	$(MPIF90) $(LDFLAGS) -o dene.x \
        $(PWOBJS) $(PHOBJS) obm_tests.o $(QEMODS) $(LIBOBJS) $(LIBS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f ../bin/turbo_lanczos.x ../bin/turbo_davidson.x *.x *.o *~ *.d *.mod *.i work.pc *.a *.F90

include make.depend

# DO NOT DELETE
