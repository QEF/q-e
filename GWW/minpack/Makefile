# Makefile for Minpack 

include ../../make.inc

MINOBJS = \
dpmpar.o \
fdjac2.o \
lmder1.o \
lmdif1.o \
mlagzo.o \
qrfac.o \
test.o \
enorm.o \
lmder.o \
lmdif.o \
lmpar.o \
mlegzo.o \
qrsolv.o

all : minpacklib.a 

minpacklib.a : $(MINOBJS)
	$(AR) $(ARFLAGS) $@ $(MINOBJS)

clean :
	- /bin/rm -f minpacklib.a *.o 

