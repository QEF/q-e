# Makefile for DAVID

include ../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG).

UTIL = clocks_handler.o \
cuda_util.o \
divide.o \
data_buffer.o \
error_handler.o \
find_free_unit.o \
fletcher32_mod.o \
mem_counter.o \
mp.o \
mp_base.o \
mp_base_gpu.o \
mp_bands_util.o \
parallel_include.o \
util_param.o \
thread_util.o

all : libutil.a


libutil.a:  $(UTIL)
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
