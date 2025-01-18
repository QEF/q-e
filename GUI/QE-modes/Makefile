TOPDIR   = $(shell cd ../..; pwd)
VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
GEN_MODE = $(TOPDIR)/dev-tools/gen-emacs-mode
LATEX    = pdflatex

ELISP_FILES = *.el *.elc

PWTK_CMDS = $(shell pwtk -s $(TOPDIR)/dev-tools/exported-cmds.pwtk)

# automatic ...
ALL_DEF_FILES = $(shell find $(TOPDIR) | grep INPUT_ | egrep def$$)

## ... manual (i.e. perhaps more secure, as to have some control)
#ALL_DEF_FILES = \
#	$(TOPDIR)/PW/Doc/INPUT_PW.def \
#	$(TOPDIR)/PWCOND/Doc/INPUT_PWCOND.def \
#	$(TOPDIR)/PP/Doc/INPUT_IMPORTEXPORT_BINARY.def \
#	$(TOPDIR)/PP/Doc/INPUT_BANDS.def \
#	$(TOPDIR)/PP/Doc/INPUT_bgw2pw.def \
#	$(TOPDIR)/PP/Doc/INPUT_PROJWFC.def \
#	$(TOPDIR)/PP/Doc/INPUT_pw2bgw.def \
#	$(TOPDIR)/PP/Doc/INPUT_PP.def \
#	$(TOPDIR)/PP/Doc/INPUT_molecularpdos.def \
#	$(TOPDIR)/PP/Doc/INPUT_DOS.def \
#	$(TOPDIR)/PHonon/Doc/INPUT_PH.def \
#	$(TOPDIR)/PHonon/Doc/INPUT_D3.def \
#	$(TOPDIR)/atomic/Doc/INPUT_LD1.def \
#	$(TOPDIR)/CPV/Doc/INPUT_CP.def \
#	$(TOPDIR)/CPV/Doc/INPUT_CPPP.def \
#	$(TOPDIR)/TDDFPT/Doc/INPUT_EELS.def \
#	$(TOPDIR)/TDDFPT/Doc/INPUT_Spectrum.def \
#	$(TOPDIR)/TDDFPT/Doc/INPUT_Lanczos.def \
#	$(TOPDIR)/TDDFPT/Doc/INPUT_Davidson.def \
#	$(TOPDIR)/West/Doc/INPUT_WFREQ.def \
#	$(TOPDIR)/West/Doc/INPUT_WSTAT.def \
#	$(TOPDIR)/NEB/Doc/INPUT_NEB.def

PDFS = user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

generate = if test ! -d qe-modes; then mkdir qe-modes; fi; cd qe-modes; $(GEN_MODE)


all: modes doc
modes: make-modes make-master
make-modes: qe pw neb ph pp ld1 cp pwtk


# generate specific modes ...

pw: qe-funcs 
	$(generate) pw

neb: qe-funcs 
	$(generate) $(TOPDIR)/NEB/Doc/INPUT_NEB.def $(TOPDIR)/PW/Doc/INPUT_PW.def

ph: qe-funcs
	$(generate) ph

pp: qe-funcs
	$(generate) pp

ld1: qe-funcs
	$(generate) -modename "QE-ld1.x (atomic)" $(TOPDIR)/atomic/Doc/INPUT_LD1.def

cp: qe-funcs
	$(generate) $(TOPDIR)/CPV/Doc/INPUT_CP.def


# generate generic QE mode

qe: 
	$(generate) -funcs -mode qe -modename QE-generic $(ALL_DEF_FILES)

pwtk: 
	$(generate) -mode pwtk -modename PWTK -template pwtk-mode.el.tcl -nmlprefix "" -cmds "$(PWTK_CMDS)" $(ALL_DEF_FILES)

# generate other required elisp files ...

qe-funcs:
	- if test ! -f qe-modes/qe-funcs.el; then $(MAKE) qe; fi

make-master:
	$(generate) -masteronly -mode qe $(ALL_DEF_FILES)

# other targets ...

doc: 
	cd Doc; $(MAKE)

clean:
	- cd Doc; $(MAKE) clean

veryclean: clean
	- rm *~ */*~ */*.bak
	- rm -rf QE-modes-$(VERSION)
	- rm *.tar
	- rm *tar.gz
	- rm -rf */auto
	- cd Doc; $(MAKE) distclean

distclean: veryclean


tar: doc
	cd Doc; $(MAKE) clean; # purify the Doc/ directory
	if test -f QE-modes.tar; then rm -f QE-modes.tar; fi
	if test -f QE-modes-$(VERSION).tar.gz; then rm -f QE-modes-$(VERSION).tar.gz; fi
	find ./ -type f | grep -v -e '\.svn' -e '/\.' -e '/auto/' -e '/qe-modes.templates' \
		-e '/Makefile' -e'~$$' -e '\.tar$$' -e '\.gz$$' -e '\.tgz$$' | xargs tar rvf QE-modes.tar
	mkdir QE-modes-$(VERSION)
	cd QE-modes-$(VERSION); tar xvf ../QE-modes.tar
	tar cvf QE-modes-$(VERSION).tar QE-modes-$(VERSION)/
	gzip QE-modes-$(VERSION).tar
	rm -rf QE-modes.tar QE-modes-$(VERSION)
