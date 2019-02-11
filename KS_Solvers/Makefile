# Makefile for  KS_Solvers
sinclude ../make.inc

# the following directive prevents execution of this makefiile in parallel
# Problem is, libks_solvers.a depends upon all libdavid libcg etc, but if
# those dependencies are added explicitly to libks_solvers.a, the library
# is relinked every time, forcing relinking of many executables; if not, 
# parallel make may execute libks_solvers.a before all other lib* complete 
.NOTPARALLEL :

default: all

all: libdavid_rci libdavid libcg libppcg libparo libnewsolver libks_solvers.a

ALLOBJS = \
Davidson/david_param.o \
Davidson/cegterg.o \
Davidson/regterg.o \
Davidson_RCI/david_rci.o \
CG/cg_param.o \
CG/ccgdiagg.o \
CG/rcgdiagg.o \
CG/rotate_wfc_gamma.o \
CG/rotate_wfc_k.o \
PPCG/ppcg_param.o \
PPCG/ppcg_gamma.o \
PPCG/ppcg_k.o
# add here other objects, e.g. ParO/*.o NewSolver/*.o

libks_solvers.a: $(ALLOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

libdavid_rci : 
	if test -d Davidson_RCI ; then \
	( cd Davidson_RCI ; $(MAKE) all || exit 1 ) ; fi

libdavid : 
	if test -d Davidson ; then \
        ( cd Davidson ; $(MAKE) all || exit 1 ) ; fi

libcg : 
	if test -d CG ; then \
        ( cd CG ; $(MAKE) all || exit 1 ) ; fi

libppcg : 
	if test -d PPCG ; then \
        ( cd PPCG ; $(MAKE) all || exit 1 ) ; fi

libparo : 
	if test -d ParO ; then \
        ( cd ParO ; $(MAKE) all || exit 1 ) ; fi

libnewsolver : 
	if test -d NewSolver ; then \
        ( cd NewSolver ; $(MAKE) all || exit 1 ) ; fi

clean : 
	if test -d Davidson_RCI ; then ( cd Davidson_RCI ; $(MAKE) clean ) ; fi
	if test -d Davidson     ; then ( cd Davidson ; $(MAKE) clean ) ; fi
	if test -d CG           ; then ( cd CG ; $(MAKE) clean ) ; fi
	if test -d PPCG         ; then ( cd PPCG ; $(MAKE) clean ) ; fi
	if test -d ParO         ; then ( cd ParO ; $(MAKE) clean ) ; fi
	if test -d NewSolver    ; then ( cd NewSolver ; $(MAKE) clean ) ; fi
	rm -f libks_solvers.a

