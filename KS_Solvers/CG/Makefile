# Makefile for DAVID

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../ELPA/src $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

CG = \
ccgdiagg.o \
rcgdiagg.o \
export_gstart_2_cg.o \
set_mpi_comm_4_cg.o \
cg_param.o \
rotate_wfc_gamma.o \
rotate_wfc_k.o 

all : libcg.a


libcg.a:  $(CG)
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
