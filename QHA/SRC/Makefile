.SUFFIXES: .f90 .o
FC = ifort
LD = $(FC) -static

Proj_x    = Atom_projected_properties.x
FQHA_x    = F_QHA.x
Ghost_x   = Ghost_DOS.x
Partial_x = Partial_phonon_DOS.x
MSD_x     = Mean_square_displacement.x
Atom_x    = atom_info.x

FFLAGS = -FR

OBJ1 =  Mean_square_displacement.o

OBJ2 =	Atom_projected_properties.o

OBJ3 =  F_QHA.o

OBJ4 =  Ghost_DOS.o

OBJ5 =  Partial_phonon_DOS.o

OBJ6 =  atom_info.o
 
all:	MSD Proj FQHA Ghost Partial Atom

MSD:	$(OBJ1)
	$(LD) -o $(MSD_x) $(OBJ1)   

Proj:	$(OBJ2)
	$(LD) -o $(Proj_x) $(OBJ2)   

FQHA:	$(OBJ3)
	$(LD) -o $(FQHA_x) $(OBJ3)  

Ghost:	$(OBJ4)
	$(LD) -o $(Ghost_x) $(OBJ4)  

Partial:$(OBJ5)
	$(LD) -o $(Partial_x) $(OBJ5)  

Atom   :$(OBJ6)
	$(LD) -o $(Atom_x) $(OBJ6)  

.f90.o : 
	$(FC) $(FFLAGS) -c  $<

clean:
	rm -f *.o         

