
CODEDIR=CPV

include	../make.sys


all: cp90

#   Name of the program
PROGNAME= cp.x

######################################################################

WRAPPERS= wrapper.o

#   Fortran objects
FOBJS= modules.o crayfft.o para.o \
read_pseudo.o cplib.o cpr.o input.o \
macdep.o which_dft.o $(WRAPPERS) \
restart.o cpr_mod.o cprsub.o fft_cp.o

######################################################################

default:
	make $(PROGNAME)

cp90 : $(FOBJS) cpflush.o
	$(LD) -o $(PROGNAME) $(FOBJS) ../Modules/*.o cpflush.o \
          ../flib/flib.a ../flib/eispack.o $(LDFLAGS)
	chmod a+r $(PROGNAME)
	chmod a+x $(PROGNAME)


clean_:
	- rm -f *.o *.i core* fort* *.mod *.s *.d work.pc cp.x

include .dependencies
include ../make.rules
