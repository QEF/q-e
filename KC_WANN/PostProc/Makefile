# Makefile for KC PP. 

include ../../make.inc

IFLAGS         = -I../../include -I../src
#LIBOBJS	       = ../../clib/clib.a ../../iotk/src/libiotk.a

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG)../../PHonon/PH \
          $(MOD_FLAG).

PP_OBJS=bcast_kc_pp_input.o \
input_pp_summary.o \
kc_pp_readin.o \
read_hr.o \
kc_bands.o 

SH_OBJS= sh_setup.o 


KC_OBJS = \
../src/kc_wann_com.o \
../src/interpolation.o \
../src/convert_kpts_names.o \
../src/self_hartree.o \
../src/setup_coulomb_exx.o \
../src/bare_pot.o \
../src/kc_io_new.o \
../src/read_wannier.o 

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
PHOBJS = ../../PHonon/PH/libph.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../upflib/libupf.a \
         ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a


TLDEPS= kc  #Dependency. To be compiled first

all : tldeps interp_r_to_k.x sh.x

interp_r_to_k.x: interp_r_to_k.o $(PP_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ interp_r_to_k.o \
        $(PP_OBJS) $(KC_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)  \
        $(LIBOBJS) $(LIBS) $(QELIBS)

sh.x: compute_self_hartree.o $(SH_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ compute_self_hartree.o \
        $(SH_OBJS) $(KC_OBJS) $(PHOBJS) $(LRMODS) $(PWOBJS) $(QEMODS)  \
        $(LIBOBJS) $(LIBS) $(QELIBS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../ ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *__genmod.f90 *__genmod.mod

clean_partial :
	- /bin/rm -f *__genmod.f90 *__genmod.mod

include make.depend
