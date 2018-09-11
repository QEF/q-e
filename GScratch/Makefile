# Makefile for GScratch

include ../make.inc

MODFLAGS=

GSCRATCH = tiny_buffer.o \
gbuffers.o 

all : libgscratch.a

libgscratch.a:  $(GSCRATCH)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@  

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

include make.depend
