.SUFFIXES: .f90 .f .o
FC = ifort
LD = $(FC) -static
Debye_x = Debye.x

FFLAGS = -O3 

SRCS  =	Debye.f90 Debye_T.f debye3.f cheval.f d1mach.f 

OBJS  =	Debye.o Debye_T.o debye3.o cheval.o d1mach.o 

Debye_x:$(OBJS)
	$(LD) -o $(Debye_x) $(OBJS)

.f.o : 
	$(FC) $(FFLAGS) -c  $<

.f90.o : 
	$(FC) $(FFLAGS)  -c  $<

clean:
	rm -f *.o *.x
