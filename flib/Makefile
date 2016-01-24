# Makefile for flib

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG).

OBJS = \
atomic_number.o \
capital.o \
cryst_to_car.o \
dost.o \
erf.o \
expint.o \
find_free_unit.o \
functionals.o \
generate_k_along_lines.o \
lsda_functionals.o \
more_functionals.o \
has_xml.o \
inpfile.o \
int_to_char.o \
latgen.o \
linpack.o \
metagga.o \
matches.o \
plot_io.o \
radial_gradients.o \
rgen.o \
recips.o \
remove_tot_torque.o \
set_hubbard_l.o \
simpsn.o \
sort.o \
sph_bes.o \
sph_dbes.o \
trimcheck.o \
test_input_file.o \
date_and_tim.o \
volume.o \
dylmr2.o \
ylmr2.o \
wgauss.o \
w0gauss.o \
w1gauss.o \
deviatoric.o


all: flib.a

flib_only: flib.a

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

clean :
	- /bin/rm -f *.a *.o *.mod *.i *.F90 core* *.L

include make.depend
