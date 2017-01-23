# Makefile for LAXlib

include ../make.inc

MODFLAGS= $(MOD_FLAG)../ELPA/src

LAX = la_types.o \
la_error.o \
dspev_drv.o \
ptoolkit.o \
transto.o \
distools.o \
zhpev_drv.o \
la_param.o


all : libqela.a


libqela.a:  $(LAX)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

la_test.x : test.o libqela.a 
	$(LD) $(LDFLAGS) -o la_test.x test.o libqela.a $(LIBS)

TEST : la_test.x

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

include make.depend
