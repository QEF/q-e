.SUFFIXES: .f90 .o

include ./make.inc

kpoints_x = k_for_bands.x
bands_x   = bands_to_gnuplot.x
Eminmax_x = E_min_max.x

OBJ1 =	K_for_bands.o

OBJ2 =	Bands_to_gnuplot.o

OBJ3 =	E_min_max.o

all:	Kpts  E_bands  Plot

Kpts:	$(OBJ1)
	$(LD) -o $(kpoints_x) $(OBJ1)  

E_bands:$(OBJ2)
	$(LD) -o $(bands_x) $(OBJ2)   

Plot:	$(OBJ3)
	$(LD) -o $(Eminmax_x) $(OBJ3)   


.f90.o : 
	$(FC) $(FFLAGS) -c  $<

clean:
	rm -f *.o         

clean_all: clean
	\rm -f *.x
