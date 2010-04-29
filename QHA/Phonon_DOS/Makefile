.SUFFIXES: .f .o
FC = ifort
LD = $(FC) -static
DOS  = phonon_dos.x 
TET  = tetra.x

#Linux
#FFLAGS = -O3 -ffast-math -fno-f2c

FFLAGS = -O3 

OBJ1 =	phonon_dos.o k_brillouin.o generate_tetra.o  det3.o det4.o \
	Tetrahedra.o Integration.o 

OBJ2 =	tetra.o k_brillouin.o generate_tetra.o  det3.o det4.o 

all:	tetra phdos

tetra:	$(OBJ2)
	$(LD) -o $(TET) $(OBJ2) $(LIBS) 

phdos:	$(OBJ1)
	$(LD) $(OBJ1) $(LIBS) -o $(DOS)

.f.o : 
	$(FC) $(FFLAGS) -c  $<

clean:
	rm -f *.o         
