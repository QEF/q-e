# Makefile for clib

include ../make.inc

OBJS = \
customize_signals.o \
stack.o \
c_mkdir.o \
copy.o \
cptimer.o \
eval_infix.o \
indici.o \
md5.o \
md5_from_file.o \
memstat.o \
ptrace.o \
sockets.o \
qmmm_aux.o \
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
