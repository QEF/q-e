# Makefile for PH
sinclude ../make.inc

# see comments in file tldeps for explanations
TLDEPS=pwlibs lrmods
include ../install/tldeps

all: ph phcg # finite_diffs

ph: phlibs
	( cd PH ; $(MAKE) all || exit 1 )

phlibs:
	( cd PH ; $(MAKE) phlibs || exit 1 )

phcg:
	( cd Gamma ; $(MAKE) all || exit 1 )

finite_diffs:
	( cd FD ; $(MAKE) all || exit 1 )

clean: ph_clean phcg_clean examples_clean finite_diffs_clean

ph_clean:
	( cd PH ; $(MAKE) clean )

phcg_clean:
	( cd Gamma ; $(MAKE) clean )

#third_order_q_clean:
#	( cd D3q ; $(MAKE) clean )

finite_diffs_clean:
	( cd FD ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
