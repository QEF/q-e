# Makefile for clib

include ../make.sys

OBJS = \
stack.o \
c_mkdir.o \
cptimer.o \
fft_stick.o \
indici.o \
memstat.o \
qsort.o

all : clib.a

clib.a : $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

source :
	co -l $(OBJS:.o=.c)

clean :
	- rm -f clib.a *.o *.mod *.i core*

include make.depend
