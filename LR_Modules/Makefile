# Makefile for LR_Modules

include ../make.sys

# location of needed modules

LRMODULES = 

TLDEPS=libfft mods

all : liblrmod.a

liblrmod.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: version all clean

include make.depend
