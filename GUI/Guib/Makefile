VERSION = $(shell cat VERSION)
TMPDIR  = /tmp/Guib-$(VERSION)
HERE    = $(shell pwd)
REPACK  = $(HERE)/sys_utils/repack.sh

# ------------------------------------------------------------------------
#  INDIVIDUAL FILES (per directory)	
# ------------------------------------------------------------------------

ROOT_FILES = \
	README THANKS NEWS INSTALL TODO VERSION ChangeLog \
	BUGS AUTHORS COPYING COPYRIGHT \
	guib  guib.bat  guib.tcl  guib.theme  init.tcl  pkgIndex.tcl 

DOC_FILES = doc/html/*.html doc/usage.pdf

EXAMPLES_FILES = \
	examples/*.tcl \
	examples/*.inp \
	examples/README \
	examples/run.sh

EXTERNAL_FILES = external/

IMAGES_FILES = images/*.gif images/icons.license

SRC_FILES = \
	src/tclIndex \
	src/aindex \
	src/*tcl 

LIB_FILES = \
	lib/tclIndex lib/pkgIndex.tcl \
	lib/tclUtils.tcl lib/tkUtils.tcl


# ------------------------------------------------------------------------
#  FILES grouped by content
# ------------------------------------------------------------------------
GUIB_FILES = \
	$(ROOT_FILES) $(EXAMPLES_FILES) $(IMAGES_FILES) $(SRC_FILES)

LIBS_FILES = $(LIB_FILES) $(EXTERNAL_FILES)

ALL_FILES = $(DOC_FILES) $(GUIB_FILES) $(LIBS_FILES) 

# make nothing by default !!!

default: 
	@echo "Try one of the following:"
	@echo ""
	@echo "make guib-archive       ... makes GUIB *.tgz archive (with documentation)"
	@echo "make guib-archive-nodoc ... makes GUIB *.tgz archive (no documentation)"
	@echo "make clean              ... rm *.tgz *.tar *.tgz *~"
	@echo "make veryclean          ... clean + rm *~ (all directories)"
	@echo ""

guib-archive:
	if test ! -d doc/usage; then cd doc; $(MAKE); fi
	tar --exclude=.svn -cvhf guib.tar $(ALL_FILES)	
	$(REPACK) guib guib-$(VERSION) Guib-$(VERSION)
	rm guib.tar

guib-archive-nodoc:
	tar --exclude=.svn -cvhf guib.tar $(GUIB_FILES) $(LIBS_FILES)
	$(REPACK) guib guib-$(VERSION) Guib-$(VERSION)
	rm guib.tar

clean:
	-rm -f *~
	-rm -f *.tar
	-rm -f *.tar.gz
	-rm -f *.tgz

veryclean: clean
	-for file in `find . | egrep ~$$`; do echo "cleaning ...  $$file"; rm -f $$file; done
