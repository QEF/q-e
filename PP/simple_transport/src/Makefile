#
#

include ../../../make.inc

# location of needed modules and included files (if any)
# MODFLAGS= $(MOD_FLAG)../../../Modules

all: ef dos fermi_int_0 fermi_int_1
clean:
	/bin/rm *.o *.x *.mod

ef: ef.o cryst_to_car.o lint.o 
	$(LD) $(LDFLAGS) cryst_to_car.o lint.o ef.o  -o ef.x

dos: dos.o smearing_mod.o cryst_to_car.o lint.o vband_ibz.o
	$(LD) $(LDFLAGS) smearing_mod.o cryst_to_car.o lint.o vband_ibz.o dos.o -o dos.x
dos.o: smearing_mod.o

fermi_int_0: fermi_int_0.o
	$(LD) $(LDFLAGS) cryst_to_car.o lint.o vband_ibz.o fermi_int_0.o -o fermi_int_0.x

fermi_int_1:  fermi_int_1.o reducegrid.o
	$(LD) $(LDFLAGS) cryst_to_car.o reducegrid.o lint.o vband_ibz.o fermi_int_1.o -o fermi_int_1.x
