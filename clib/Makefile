
include ../make.sys

all: clib.a

OBJS = c_getenv.o c_today.o cpflush.o factor.o indici.o \
       memstat.o qsort.o readocc.o round2.o cptimer.o fft_stick.o 


clib.a : $(OBJS) 
	ar $(ARFLAGS) $@ $?
	$(RANLIB) $@

source:
	co -l $(OBJS:.o=.c)

clean_:
	- rm -f clib.a *.o *.mod *.i core*

include ../make.rules
