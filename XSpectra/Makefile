# Makefile for X-ray spectra calculation (XSPECTRA)

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../EE $(MOD_FLAG)../PW $(MOD_FLAG)../GIPAW $(MOD_FLAG).

XOBJS = \
./xspectra_mod.o \
./radin_mod.o \
./mygetK.o \
./read_file_xspectra.o \
./ipoolscatter.o \
./init_xspectra_ldau.o \
./lr_sm1_psi.o

PPOBJS = \
../PP/openfil_pp.o\
../PP/stop_pp.o

GIPAWOBJS=../GIPAW/paw_gipaw.o \
../GIPAW/init_gipaw_1.o \
../GIPAW/init_gipaw_2.o

QEMODS = ../Modules/libqemod.a
PWOBJS = ../PW/libpw.a 
EEOBJS = ../EE/libee.a

all : xspectra.x 

xspectra.x :   xspectra.o $(XOBJS) $(GIPAWOBJS) $(QEMODS) $(PWOBJS)
	$(MPIF90) -o $@  xspectra.o $(XOBJS) $(PPOBJS) $(GIPAWOBJS) \
		$(PWOBJS) $(EEOBJS) $(QEMODS) $(LDFLAGS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../XSpectra/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L

include make.depend
# DO NOT DELETE
