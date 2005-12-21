# Makefile for CPVIB/PWVIB

.PHONY : cpvib.x pwvib.x

all : cpvib pwvib

cpvib   : cpvib.x
cpvib.x : 
	$(MAKE) -f makefile.cpvib

pwvib   : pwvib.x
pwvib.x : 
	$(MAKE) -f makefile.pwvib

clean :
	- /bin/rm -f cpvib.x pwvib.x *.o *.mod version.h *.i core* *.F90 fort* \
		*.cpp *.d work.pc  *.s

