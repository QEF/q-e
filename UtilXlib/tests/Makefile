# Makefile for DAVID

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG). $(MOD_FLAG)../

DFLAGS=$$DFLAGMATRIX

SRCS = test_mp_count_nodes.f90 \
       test_mp_bcast_i1.f90 \
       test_mp_bcast_i1_gpu.f90 \
       test_mp_bcast_iv.f90 \
       test_mp_bcast_iv_gpu.f90 \
       test_mp_bcast_im.f90 \
       test_mp_bcast_im_gpu.f90 \
       test_mp_bcast_it.f90 \
       test_mp_bcast_it_gpu.f90 \
        autotest_mp_bcast_i1.f90\
        autotest_mp_bcast_iv.f90\
        autotest_mp_bcast_im.f90\
        autotest_mp_bcast_it.f90\
        autotest_mp_bcast_r1.f90\
        autotest_mp_bcast_rv.f90\
        autotest_mp_bcast_rm.f90\
        autotest_mp_bcast_rt.f90\
        autotest_mp_bcast_r1.f90\
        autotest_mp_bcast_rv.f90\
        autotest_mp_bcast_rm.f90\
        autotest_mp_bcast_rt.f90\
        autotest_mp_bcast_c1.f90\
        autotest_mp_bcast_cv.f90\
        autotest_mp_bcast_cm.f90\
        autotest_mp_bcast_ct.f90\
        autotest_mp_bcast_c1.f90\
        autotest_mp_bcast_cv.f90\
        autotest_mp_bcast_cm.f90\
        autotest_mp_bcast_ct.f90\
        autotest_mp_bcast_gpu_i1.f90\
        autotest_mp_bcast_gpu_iv.f90\
        autotest_mp_bcast_gpu_im.f90\
        autotest_mp_bcast_gpu_it.f90\
        autotest_mp_bcast_gpu_r1.f90\
        autotest_mp_bcast_gpu_rv.f90\
        autotest_mp_bcast_gpu_rm.f90\
        autotest_mp_bcast_gpu_rt.f90\
        autotest_mp_bcast_gpu_r1.f90\
        autotest_mp_bcast_gpu_rv.f90\
        autotest_mp_bcast_gpu_rm.f90\
        autotest_mp_bcast_gpu_rt.f90\
        autotest_mp_bcast_gpu_c1.f90\
        autotest_mp_bcast_gpu_cv.f90\
        autotest_mp_bcast_gpu_cm.f90\
        autotest_mp_bcast_gpu_ct.f90\
        autotest_mp_bcast_gpu_c1.f90\
        autotest_mp_bcast_gpu_cv.f90\
        autotest_mp_bcast_gpu_cm.f90\
        autotest_mp_bcast_gpu_ct.f90\

EXECS = $(SRCS:.f90=.x)

all: common $(EXECS)

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

