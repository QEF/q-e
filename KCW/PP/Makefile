# Makefile for KC PP. 

include ../../make.inc

VPATH = @srcdir@

IFLAGS         = -I../../include -I../src

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../LR_Modules \
          $(MOD_FLAG).

PP_OBJS=bcast_kcw_pp_input.o \
input_pp_summary.o \
kcw_pp_readin.o \
read_hr.o \
kcw_bands.o\

SH_OBJS= sh_setup.o\
input_pp_summary.o 

KC_OBJS = ../src/libkcw.a

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../upflib/libupf.a ../../XClib/xc_lib.a \
         ../../FFTXlib/src/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a \
          ../..//MBD/libmbd.a

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS= pwlibs lrmods kcwlib
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../../install/tldeps

all : kcwpp_sh.x kcwpp_interp.x clean_partial merge_wann.x

merge_wann.x: merge_wann.o
	$(LD) $(LDFLAGS) -o $@ merge_wann.o 
	- ( cd ../../bin ; ln -fs ../KCW/PP/$@ . )

kcwpp_interp.x: interp_r_to_k.o $(PP_OBJS) $(LRMODS) $(PWOBJS) $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ interp_r_to_k.o \
        $(PP_OBJS) $(KC_OBJS) $(LRMODS) $(PWOBJS) $(QEMODS)  \
        $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../KCW/PP/$@ . )

kcwpp_sh.x:  compute_self_hartree.o $(SH_OBJS) $(LRMODS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ compute_self_hartree.o \
        $(SH_OBJS) $(KC_OBJS) $(LRMODS) $(PWOBJS) $(QEMODS)  \
        $(LIBS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../KCW/PP/$@ . )

ups.x: ups.o $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ ups.o $(PWOBJS) $(QEMODS) \
	$(LIBS) $(QELIBS)

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *__genmod.f90 *__genmod.mod

clean_partial :
	- /bin/rm -f *__genmod.f90 *__genmod.mod


include make.depend


