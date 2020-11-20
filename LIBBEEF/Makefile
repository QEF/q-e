# Makefile for libbeef

include ../make.inc

OBJS = \
pbecor.o \
beefun.o

# PGI requires -c11 to understand __thread
# See https://www.pgroup.com/userforum/viewtopic.php?t=5504
ifneq (,$(findstring -D__PGI,$(DFLAGS)))
	CFLAGS := -c11 $(CFLAGS)
endif

all : libbeef.a

libbeef.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

source :
	co -l $(OBJS:.o=.c)

clean :
	- rm -f libbeef.a *.o *.mod *.i core*

include make.depend
