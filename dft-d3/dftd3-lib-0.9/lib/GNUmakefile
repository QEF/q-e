################################################################################
#
# Library makefile
#
# Needs following variables
#
# FC  Fortran compiler
# FCFLAGS  Fortran compiler flags
# LN  Linker
# LNFLAGS  Linker flags
# SRCDIR  Folder, where source files can be found
#
################################################################################
.SUFFIX:
.SUFFIX: .f90 .o

LIBD3 = libdftd3.a

all: $(LIBD3)

vpath % $(SRCDIR)

.PHONY: clean distclean


OBJS  = api.o common.o core.o pars.o sizes.o

$(LIBD3): $(OBJS)
	ar r $@ $^

clean:
	rm -f *.o *.mod *.a

distclean: clean
	rm -f *.mod $(LIBD3)


%.o: %.f90
	$(FC) $(FCFLAGS) -c $< -o $@


# Dependencies
api.o: common.o core.o sizes.o
common.o:
core.o: common.o pars.o sizes.o
pars.o: common.o sizes.o
sizes.o:
