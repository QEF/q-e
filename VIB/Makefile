# Makefile for CPVIB/PWVIB

.PHONY : cpvib.x pwvib.x

include ../make.sys

TLDEPS=bindir mods libs libiotk pw cp

all : tldeps cpvib.x pwvib.x

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :
 
cpvib.x : 
	$(MAKE) -f makefile.cpvib DFLAGS='$(DFLAGS) -DDFT_CP'

pwvib.x : 
	$(MAKE) -f makefile.pwvib DFLAGS='$(DFLAGS) -DDFT_PW' 

clean :
	- /bin/rm -f cpvib.x pwvib.x *.o *.mod version.h *.i core* *.F90 fort* \
		*.cpp *.d work.pc  *.s

