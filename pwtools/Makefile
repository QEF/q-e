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
         ../PW/irrek.o ../PW/mode_group.o \
         ../PW/cft_3.o ../PW/error_handler.o
MODULES = ../Modules/parameters.o
#
# Targets
#
all: band_plot.x q2r.x matdyn.x dynmat.x fqha.x

band_plot.x:  band_plot.o
	$(LD) -o band_plot.x band_plot.o  $(LFLAGS)

q2r.x: q2r.o rigid.o
	$(LD) -o q2r.x q2r.o rigid.o $(PWOBJS) $(MODULES) $(LFLAGS)

matdyn.x: matdyn.o rigid.o
	$(LD) -o matdyn.x matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LFLAGS)

dynmat.x: dynmat.o rigid.o
	$(LD) -o dynmat.x dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LFLAGS)

fqha.x:  fqha.o
	$(LD) -o fqha.x fqha.o  $(LFLAGS)

clean_:
	rm -f *.x *.o *~ *.F90 *.mod *.d *.i work.pc
