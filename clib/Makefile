# Makefile for clib

include ../make.sys
include ../make.rules

OBJS = \
c_mkdir.o \
cpflush.o \
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
