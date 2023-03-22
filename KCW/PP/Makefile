# Makefile for KC PP. 

include ../../make.inc

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
kcw_bands.o 

SH_OBJS= sh_setup.o 

KC_OBJS = ../src/libkcw.a

LRMODS = ../../LR_Modules/liblrmod.a
PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a ../../upflib/libupf.a ../../XClib/xc_lib.a \
         ../../FFTXlib/src/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a \
          ../..//MBD/libmbd.a

TLDEPS= kcw-lib  #Dependency. To be compiled first

all : tldeps kcwpp_sh.x kcwpp_interp.x clean_partial

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
	$(LD) $(LDFLAGS) -o $@ ups.o ../../PP/src/stop_pp.o $(PWOBJS) $(QEMODS) \
	$(LIBS) $(QELIBS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../ ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L *__genmod.f90 *__genmod.mod

clean_partial :
	- /bin/rm -f *__genmod.f90 *__genmod.mod

include make.depend
