#
# Makefile for Modules
#
MODULES= kind.o \
parallel_include.o \
shmem_include.o \
parameters.o \
io_base.o  \
io_global.o  \
mp.o  \
mp_global.o  \
mp_wave.o  \
parser.o \
readpseudo.o  \
pseudodata.o  \
pseudo_types.o
#
include ../make.rules
#
include ../make.sys
#
# targets
#
all: $(MODULES)
#
clean:
	-/bin/rm -f *.o *.d *~ *.F90 *.mod
veryclean: clean
	-/bin/rm .dependencies
#
include .dependencies

