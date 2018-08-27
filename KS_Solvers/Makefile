# Makefile for  KS_Solvers
sinclude ../make.inc

default: all

all: libks_solvers.a

libks_solvers.a: libdavid_rci libdavid libcg libppcg libparo libnewsolver
	$(AR) $(ARFLAGS) $@ */*.o       
	$(RANLIB) $@    

libdavid_rci : 
	if test -d Davidson_RCI ; then \
	( cd Davidson_RCI ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

libdavid : 
	if test -d Davidson ; then \
        ( cd Davidson ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

libcg : 
	if test -d CG ; then \
        ( cd CG ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

libppcg : 
	if test -d PPCG ; then \
        ( cd PPCG ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

libparo : 
	if test -d ParO ; then \
        ( cd ParO ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

libnewsolver : 
	if test -d NewSolver ; then \
        ( cd NewSolver ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

clean : 
	if test -d Davidson_RCI ; then ( cd Davidson_RCI ; $(MAKE) clean ) ; fi
	if test -d Davidson     ; then ( cd Davidson ; $(MAKE) clean ) ; fi
	if test -d CG           ; then ( cd CG ; $(MAKE) clean ) ; fi
	if test -d PPCG         ; then ( cd PPCG ; $(MAKE) clean ) ; fi
	if test -d ParO         ; then ( cd ParO ; $(MAKE) clean ) ; fi
	if test -d NewSolver    ; then ( cd NewSolver ; $(MAKE) clean ) ; fi
	rm -f libks_solvers.a

