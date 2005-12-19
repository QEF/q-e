# Makefile for CPVIB/PWVIB

all : cpvib.x pwvib.x

cpvib   : cpvib.x
cpvib.x : 
	$(MAKE) -f makefile.cpvib

pwvib   : pwvib.x
pwvib.x : 
	$(MAKE) -f makefile.pwvib

clean :
	- /bin/rm -f cpvib.x pwvib.x *.o *.mod version.h *.i core* *.F90 fort* \
		*.cpp *.d work.pc  *.s

