#
# Makefile for tools
#
include ../make.sys
include ../make.rules
#
PWOBJS = ../PW/error.o ../PW/latgen.o ../PW/volume.o ../PW/recips.o \
         ../PW/cryst_to_car.o ../PW/cubicsym.o ../PW/hexsym.o \
         ../PW/kpoint_grid.o ../PW/sgama.o ../PW/sgam_at.o  \
         ../PW/sgam_ph.o ../PW/coset.o ../PW/multable.o \
         ../PW/smallg_q.o ../PW/dsum.o ../PW/trnvecc.o ../PW/invmat.o \
         ../PW/checksym.o ../PW/eqvect.o ../PW/setv.o  \
         ../PW/irrek.o ../PW/mode_group.o ../PW/hpsort.o ../PW/invmat.o \
         ../PW/cft_3.o ../PW/allocate.o ../PW/error_handler.o
#
# Targets
#
all: band_plot.x q2r.x matdyn.x dynmat.x band_plot.x fqha.x

band_plot.x:  band_plot.o
	$(LD) -o band_plot.x band_plot.o  $(LFLAGS)

q2r.x: q2r.o rigid.o
	$(LD) -o q2r.x q2r.o rigid.o $(PWOBJS) $(LFLAGS)

matdyn.x: matdyn.o rigid.o
	$(LD) -o matdyn.x matdyn.o rigid.o $(PWOBJS) $(LFLAGS)

dynmat.x: dynmat.o rigid.o
	$(LD) -o dynmat.x dynmat.o rigid.o $(PWOBJS) $(LFLAGS)

dos.x:  dos.o
	$(LD) -o dos.x dos.o  $(LFLAGS)

fqha.x:  fqha.o
	$(LD) -o fqha.x fqha.o  $(LFLAGS)

clean:
	rm -f *.x *.o  *.mod *.d *.F90
