# Makefile for XClib

include ../make.inc

MODFLAGS= $(MOD_FLAG)$(TOPDIR)/upflib/

XCL = \
xclib_utils_and_para.o \
xclib_error.o \
qe_kind.o \
qe_constants.o \
dft_setting_params.o \
dft_setting_routines.o \
qe_dft_list.o \
qe_funct_corr_lda_lsda.o \
qe_funct_exch_lda_lsda.o \
qe_funct_corr_gga.o \
qe_funct_exch_gga.o \
qe_funct_mgga.o \
xc_beef_interface.o \
qe_drivers_lda_lsda.o \
qe_drivers_d_lda_lsda.o \
qe_drivers_gga.o \
qe_drivers_d_gga.o \
qe_drivers_mgga.o \
xc_wrapper_d_gga.o \
xc_wrapper_d_lda_lsda.o \
xc_wrapper_gga.o \
xc_wrapper_lda_lsda.o \
xc_wrapper_mgga.o \
xc_lib.o

BEEF = \
pbecor.o \
beefun.o

XCT = \
qe_dft_refs.o

TLDEPS = libupf

all: xc_lib.a

xc_lib.a: $(XCL) $(BEEF)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

xc_infos.x : xc_infos.o qe_dft_refs.o qe_dft_list.o xc_lib.a
	$(LD) $(LD_LIBS) $(LIBXC_LIBS) $(LDFLAGS) $(DFLAGS) -o $@ \
		xc_infos.o qe_dft_refs.o qe_dft_list.o xc_lib.a $(LD_LIBS) $(LIBXC_LIBS)
	- ( cd ../bin ; ln -fs ../XClib/$@ . )

xclib_test.x : xclib_test.o $(XCT) xc_lib.a
	$(LD) ../upflib/libupf.a $(LD_LIBS) $(LIBXC_LIBS) $(LDFLAGS)  $(DFLAGS) -o $@ \
		xclib_test.o $(XCT) xc_lib.a ../upflib/libupf.a $(BLAS_LIBS) $(LD_LIBS) $(LIBXC_LIBS)
	- ( cd ../bin ; ln -fs ../XClib/$@ . )

infos : xc_infos.x

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd .. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

test : tldeps xclib_test.x

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend

# DO NOT DELETE
