# Makefile for clib

include ../make.sys
include ../make.rules

OBJS = \
c_mkdir.o \
c_today.o \
cpflush.o \
cptimer.o \
factor.o \
fft_stick.o \
indici.o \
memstat.o \
qsort.o \
readocc.o \
round2.o

all : clib.a

clib.a : $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?

source :
	co -l $(OBJS:.o=.c)

clean :
	- rm -f clib.a *.o *.mod *.i core*
