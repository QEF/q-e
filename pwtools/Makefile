#
# Makefile for tools
#
include ../make.sys
include ../make.rules
#
PWOBJS = ../PW/random.o ../PW/error.o ../PW/wsweight.o\
         ../PW/cryst_to_car.o ../PW/cubicsym.o ../PW/hexsym.o \
         ../PW/kpoint_grid.o ../PW/sgama.o ../PW/sgam_at.o  \
         ../PW/sgam_ph.o ../PW/coset.o ../PW/multable.o \
         ../PW/smallg_q.o ../PW/dsum.o ../PW/trnvecc.o ../PW/invmat.o \
         ../PW/checksym.o ../PW/eqvect.o  ../PW/cryst_to_car.o \
         ../PW/irrek.o ../PW/mode_group.o \
         ../PW/cft_3.o ../PW/error_handler.o
MODULES = ../Modules/parameters.o ../Modules/kind.o ../Modules/fft_scalar.o \
          ../Modules/constants.o

#
# Targets
#
all: band_plot.x q2r.x matdyn.x dynmat.x fqha.x dist.x ev.x kpoints.x \
     pwi2xsf pwi2xsf_new

band_plot.x:  band_plot.o
	$(LD) -o band_plot.x band_plot.o  $(LDFLAGS)

dist.x:  dist.o
	$(LD) -o dist.x dist.o $(PWOBJS) $(MODULES) $(LDFLAGS)

q2r.x: q2r.o rigid.o
	$(LD) -o q2r.x q2r.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)

matdyn.x: matdyn.o rigid.o
	$(LD) -o matdyn.x matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)

dynmat.x: dynmat.o rigid.o
	$(LD) -o dynmat.x dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LDFLAGS)

fqha.x:  fqha.o
	$(LD) -o fqha.x fqha.o  $(LDFLAGS)

ev.x:	ev.o
	$(LD) -o ev.x ev.o $(PWOBJS) $(MODULES) $(LDFLAGS)

kpoints.x: kpoints.o
	$(LD) -o kpoints.x kpoints.o $(PWOBJS) $(MODULES) $(LDFLAGS)

pwi2xsf: pwi2xsf.o pwi2xsflib.o
	$(LD) -o pwi2xsf pwi2xsf.o pwi2xsflib.o $(PWOBJS) $(MODULES) $(LDFLAGS)

pwi2xsf_new: pwi2xsf_new.o pwi2xsflib.o
	$(LD) -o pwi2xsf_new pwi2xsf_new.o pwi2xsflib.o \
        $(PWOBJS) $(MODULES) $(LDFLAGS)

clean_:
	rm -f pwi2xsf pwi2xsf_new *.x *.o *~ *.F90 *.mod *.d *.i work.pc
