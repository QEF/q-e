include ../make.sys 
 
OBJS = \
./pded.o      \
./secd.o        \
./buildAd.o     \
./buildBd.o     \
./buildGd.o     \
./buildPd.o     \
./cgd.o     \
./cgmgd.o     \
./cgmgdrvd.o     \
./maind.o     \
./mainsubd.o     \
./matvecd.o     \
./mgcsd.o     \
./mgdrvd.o    \
./mgfasd.o     \
./mgsubd.o     \
./mikpckd.o     \
./mlinpckd.o     \
./ncgd.o     \
./ninterpd.o     \
./nsmoothd.o     \
./nsord.o     \
./nwjd.o     \
./smoothd.o     \
./wjd.o     \

all : mglib.a

mglib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@


clean :
	- /bin/rm -f *.a *.o *~ *.F90 *.d *.mod *.i work.pc 
