# Makefile for PH
# Adapted from TDDFPT main Makefile

default: all

all: phonon phgamma_only third_order

phonon:
	if test -d PH ; then \
	( cd PH ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi
phgamma_only:
	if test -d Gamma; then \
	( cd Gamma ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi
third_order:
	if test -d D3; then \
	( cd D3 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

clean: phonon_clean phgamma_only_clean third_order_clean

phonon_clean:
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi

phgamma_only_clean:
	if test -d Gamma ; then \
	( cd Gamma ; if test "$(MAKE)" = ""; then make clean; \
	else $(MAKE) clean ; fi ) ; fi

third_order_clean:
	if test -d D3 ; then \
	( cd D3 ; if test "$(MAKE)" = ""; then make clean; \
	else $(MAKE) clean ; fi ) ; fi

doc:
	if test -d Doc ; then \
	(cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all ; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi
doc_clean:
	if test -d Doc ; then \
	(cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) clean ; \
	else $(MAKE) $(MFLAGS) clean ; fi ) ; fi

distclean: clean doc_clean
