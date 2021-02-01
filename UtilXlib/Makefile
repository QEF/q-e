# Makefile for DAVID

include ../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG).

UTIL = clocks_handler.o \
device_helper.o \
divide.o \
data_buffer.o \
error_handler.o \
export_gstart_2_solvers.o \
find_free_unit.o \
fletcher32_mod.o \
mem_counter.o \
mp.o \
mp_base.o \
mp_base_gpu.o \
mp_bands_util.o \
parallel_include.o \
set_mpi_comm_4_solvers.o \
util_param.o \
thread_util.o \
nvtx_wrapper.o 

all : libutil.a


libutil.a:  $(UTIL)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
