# Makefile for DAVID

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG). $(MOD_FLAG)../

DFLAGS=$$DFLAGMATRIX
LDFLAGS:=$(LDFLAGS) $(DFLAGS)

SRCS = test_mp_count_nodes.f90 \
       test_mp_bcast_i1.f90 \
       test_mp_bcast_i1_gpu.f90 \
       test_mp_bcast_iv.f90 \
       test_mp_bcast_iv_gpu.f90 \
       test_mp_bcast_im.f90 \
       test_mp_bcast_im_gpu.f90 \
       test_mp_bcast_it.f90 \
       test_mp_bcast_it_gpu.f90 \
       test_mp_bcast_iv_buffer.f90 \
       test_mp_bcast_iv_buffer_gpu.f90 \
       test_mp_bcast_lv_buffer.f90 \
       test_mp_bcast_lv_buffer_gpu.f90 \
       test_mp_bcast_rv_buffer.f90 \
       test_mp_bcast_rv_buffer_gpu.f90 \
       test_mp_max_iv_buffer_gpu.f90 \
       test_mp_max_iv_buffer.f90 \
       test_mp_max_rv_buffer.f90 \
       test_mp_max_rv_buffer_gpu.f90 \
       test_mp_min_iv_buffer_gpu.f90 \
       test_mp_min_iv_buffer.f90 \
       test_mp_min_rv_buffer.f90 \
       test_mp_min_rv_buffer_gpu.f90 \
       test_mp_sum_iv_buffer.f90 \
       test_mp_sum_rv_buffer.f90 \
       test_mp_sum_iv_buffer_gpu.f90 \
       test_mp_sum_rv_buffer_gpu.f90

# include auto generate source files. This will only work with GNU Make!
sinclude ./autotest.inc

EXECS = $(SRCS:.f90=.x)

all: common $(EXECS)

generate:
	python gen_tests.py

common: tester.o mp_world.o utils.o

%.x: %.o
	$(LD) $(LDFLAGS) $< tester.o mp_world.o utils.o -o $@ ../libutil.a


clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

