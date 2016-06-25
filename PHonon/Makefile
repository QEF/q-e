# Makefile for PH
sinclude ../make.inc

default: all

#all: phonon phgamma_only third_order third_order_q
all: phonon phgamma_only third_order finite_diffs

phonon:
	( cd PH ; $(MAKE) all || exit 1 )

phgamma_only:
	( cd Gamma ; $(MAKE) all || exit 1 )

third_order: phonon
	( cd D3 ; $(MAKE) all || exit 1 )

finite_diffs:
	( cd FD ; $(MAKE) all || exit 1 )

#third_order_q:
#	( cd D3q ; $(MAKE) all || exit 1 )

clean: phonon_clean phgamma_only_clean third_order_clean examples_clean finite_diffs_clean #third_order_q_clean

phonon_clean:
	( cd PH ; $(MAKE) clean )

phgamma_only_clean:
	( cd Gamma ; $(MAKE) clean )

third_order_clean:
	( cd D3 ; $(MAKE) clean )

#third_order_q_clean:
#	( cd D3q ; $(MAKE) clean )

finite_diffs_clean:
	( cd FD ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
