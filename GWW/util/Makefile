#Makefile for various utilities

include ../../make.inc


all: graph.x abcoeff_to_eps.x memory_pw4gww.x

graph.x: graph.f90
	$(FC) -o graph.x graph.f90
	(cd ../../bin; ln -f ../GWW/util/graph.x)

abcoeff_to_eps.x: abcoeff_to_eps.f90
	$(FC) -o abcoeff_to_eps.x abcoeff_to_eps.f90
	(cd ../../bin; ln -f ../GWW/util/abcoeff_to_eps.x)
memory_pw4gww.x: memory_pw4gww.f90
	$(FC) -o memory_pw4gww.x memory_pw4gww.f90
	(cd ../../bin; ln -f ../GWW/util/memory_pw4gww.x)

clean:
	- /bin/rm -f *.x


