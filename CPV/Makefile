
CODEDIR=CPV

include	../make.sys


all: cp90

#   Name of the program
PROGNAME= cp.x

######################################################################

WRAPPERS= wrapper.o

#   Fortran objects
FOBJS= modules.o $(FFT) para.o \
read_pseudo.o cplib.o cpr.o input.o sort.o \
macdep.o which_dft.o $(WRAPPERS) \
restart.o cpr_mod.o cell_module.o cprsub.o cpflush.o

######################################################################

default:
	make $(PROGNAME)

cp90 : $(FOBJS)
	$(LNK) -o $(PROGNAME) $(LNKFLAGS) $(FOBJS) ../Modules/*.o cpflush.o \
          ../flib/eispack.o $(LIB)
	chmod a+r $(PROGNAME)
	chmod a+x $(PROGNAME)


clean_:
	- rm -f *.o *.i core* fort* *.mod *.s *.d work.pc

include .dependencies
include ../make.rules
