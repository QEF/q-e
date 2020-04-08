#Makefile for various utilities

include ../../make.inc


all: graph.x abcoeff_to_eps.x memory_pw4gww.x

graph.x: graph.f90
	$(FC) -o graph.x graph.f90

abcoeff_to_eps.x: abcoeff_to_eps.f90
	$(FC) -o abcoeff_to_eps.x abcoeff_to_eps.f90

memory_pw4gww.x: memory_pw4gww.f90
	$(FC) -o memory_pw4gww.x memory_pw4gww.f90

clean:
	- /bin/rm -fv *.x

