# what is our Tcl/Tk version (see Make.include) ?
include Make.include

HERE    = $(shell pwd)
REPACK  = $(HERE)/sys_utils/repack.sh
STARKIT = $(HERE)/sys_utils/starkit.sh
SCRATCH = .scratch

PWGUI_VERSION    = $(shell cat VERSION)
PWGUI_TARNAME    = PWgui-$(PWGUI_VERSION)
PWGUI_DIRNAME    = PWgui-$(PWGUI_VERSION)
PWGUI_VFS        = pwgui_vfs

ROOT_FILES = \
	pwgui pwgui.bat pwgui.tcl init.tcl pwgui.settings \
	pwgui_reformat pwgui_reformat.tcl \
	README README.developer INSTALL NEWS BUGS THANKS AUTHORS \
	ChangeLog TODO COPYRIGHT COPYING VERSION
#Makefile header.tcl

#
# per directory files ... (alphabetic order)
#
DOC_FILES = \
	doc/pwdocs/INPUT*.html \
	doc/pwdocs/user_guide/

EXAMPLES_FILES = \
	examples/pw/*.pw.inp examples/pw/*.in \
	examples/neb.dat/*.dat \
	examples/ph/*.ph.inp \
	examples/pp/*.pp.inp \
	examples/projwfc/*.pr.inp \
	examples/d3/*.d3.inp \
	examples/atomic/*.in	

IMAGES_FILES = \
	images/*.gif \
	images/license.terms

MODULES_FILES = \
	modules/*/*.tcl

SRC_FILES = \
	src/pwscf.theme \
	src/about.html \
	src/aindex \
	src/tclIndex \
	src/*.itcl

TEST_FILES = tests/run-tests.sh

PWGUI_SRC_FILES = \
	$(ROOT_FILES) $(DOC_FILES) $(EXAMPLES_FILES) \
	$(IMAGES_FILES) $(MODULES_FILES) $(SRC_FILES) $(TEST_FILES)

_GUIB_DIR      = $(shell cd ../Guib; pwd)
_GUIB_VERSION  = $(shell GUIB=$(_GUIB_DIR); export GUIB; $(_GUIB_DIR)/guib -v)
_GUIB_TGZNAME  = guib-$(_GUIB_VERSION).tgz
_GUIB_DIRNAME  = Guib-$(_GUIB_VERSION)

# nothing by default

default: 
	@cat make.usage

clean:
	-rm -f *~ */*~ */*/*~
	-rm -f *.tgz 
	-rm -f *.tar

veryclean: 
	$(MAKE) clean
	-cd doc/pwdocs; $(MAKE) clean
	-cd $(PWGUI_VFS); $(MAKE) clean
	-for file in `find . | egrep ~$$`; do echo "cleaning ...  $$file"; rm -f $$file; done 

distclean: veryclean
	-cd lib/; rm -f *
	cd $(PWGUI_VFS); $(MAKE) distclean	

svninit:
	-cd doc/pwdocs; $(MAKE)
	-if test ! -d lib; then mkdir lib; fi
	-cd lib/; rm -f Guib-*;
	-cd lib/; ln -sf $(_GUIB_DIR) $(_GUIB_DIRNAME)

alldist: 
	$(MAKE) pwgui-source
	$(MAKE) pwgui-starkit-vfs

pwgui-source:  \
	_message-source \
	_create_pwgui_tar \
	_add_guib \
	_repack_pwgui

pwgui-starkit-vfs: \
	_message-starkit_vfs \
	_create_pwgui_tar \
	_add_guib \
	_create_starkit_vfs


# -- internal tagrets --

_create_pwgui_tar:
	cd doc/pwdocs; $(MAKE)
	tar cvf  pwgui.tar $(PWGUI_SRC_FILES)

_add_guib: _clean_scratch
	cd $(_GUIB_DIR); $(MAKE) guib-archive-nodoc
	\
	mv lib $(SCRATCH)
	(mkdir lib; cd lib; tar zxvf $(_GUIB_DIR)/$(_GUIB_TGZNAME))
	tar rvfh pwgui.tar lib/$(_GUIB_DIRNAME)
	rm -rf lib/
	mv $(SCRATCH) lib
	cd $(_GUIB_DIR); $(MAKE) clean

_repack_pwgui:
	$(REPACK) pwgui $(PWGUI_TARNAME) $(PWGUI_DIRNAME)
	rm pwgui.tar

_create_starkit_vfs: _clean_scratch
	cd $(PWGUI_VFS); touch make.inc make.versions; make distclean
	TOPDIR=$(HERE) PWGUI_VFS=$(PWGUI_VFS) \
	PWGUI_VERSION=$(PWGUI_VERSION) PWGUI_DIRNAME=$(PWGUI_DIRNAME) \
	GUIB_VERSION=$(_GUIB_VERSION) GUIB_DIRNAME=$(_GUIB_DIRNAME) $(STARKIT)

_clean_scratch:
	if test -d $(SCRATCH) ; then rm -rf $(SCRATCH); fi

_message-source: _check_guib
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui-source package"
	@echo "#"
	@echo "# ------------------------------------------------------------"

_message-starkit_vfs: _check_guib
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui strakit VFS (StandAlone Runtime tclKIT)      "
	@echo "#                             =    =     =          ===       "
	@echo "# ------------------------------------------------------------"

_check_guib:
	@if test -z "$(_GUIB_VERSION)" ; then \
		echo ""; echo "WARNING: GUIB package not found ..."; echo ""; exit 1; \
	fi	
