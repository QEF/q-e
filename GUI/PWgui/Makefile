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
	README README.developer INSTALL NEWS BUGS THANKS AUTHORS \
	ChangeLog TODO COPYRIGHT COPYING VERSION

#
# per directory files ... (alphabetic order)
#
DOC_FILES = \
	doc/INPUT*.html \
	doc/user_guide/

EXAMPLES_FILES = examples/

IMAGES_FILES = \
	images/*.gif \
	images/license.terms

MODULES_FILES = modules/

SRC_FILES = src/

TEST_FILES = tests/

PWGUI_SRC_FILES = \
	$(ROOT_FILES) $(DOC_FILES) $(EXAMPLES_FILES) \
	$(IMAGES_FILES) $(MODULES_FILES) $(SRC_FILES) $(TEST_FILES)

_GUIB_DIR      = $(shell cd ../Guib; pwd)
_GUIB_VERSION  = $(shell GUIB=$(_GUIB_DIR); export GUIB; $(_GUIB_DIR)/guib -v)
_GUIB_TGZNAME  = guib-$(_GUIB_VERSION).tgz
_GUIB_DIRNAME  = Guib-$(_GUIB_VERSION)


# default target == print make usage

default: 
	@cat make.usage

clean:
	-rm -f *~ */*~ */*/*~
	-rm -f *.tgz 
	-rm -f *.tar
	-cd doc; $(MAKE) clean

veryclean: clean
	-cd $(PWGUI_VFS); $(MAKE) clean
	-for file in `find . | egrep ~$$`; do echo "cleaning ...  $$file"; rm -f $$file; done 

distclean: veryclean
	-cd lib/; rm -f *
	cd $(PWGUI_VFS); $(MAKE) distclean	

init:
	-cd doc/; $(MAKE)
	-if test ! -d lib; then mkdir lib; fi
	-cd lib/; rm -f Guib-*; ln -sf $(_GUIB_DIR) $(_GUIB_DIRNAME)
	-( cd ../../bin && ln -fs ../GUI/PWgui/pwgui . )

alldist: pwgui-source pwgui-starkit-vfs

pwgui-source: init \
	_check_guib \
	_message-source \
	_create_pwgui_tar \
	_add_guib \
	_repack_pwgui

pwgui-starkit-vfs: init \
	_check_guib \
	_message-starkit_vfs \
	_create_pwgui_tar \
	_add_guib \
	_create_starkit_vfs


# -- internal tagrets --

_create_pwgui_tar: init
	cd doc/; $(MAKE) all
	tar cvf pwgui.tar $(PWGUI_SRC_FILES)

_add_guib: _clean_scratch
	cd $(_GUIB_DIR); $(MAKE) guib-archive-nodoc
	\
	mv lib $(SCRATCH)
	(mkdir lib; cd lib; tar zxvf $(_GUIB_DIR)/$(_GUIB_TGZNAME))
	tar rvfh pwgui.tar lib/$(_GUIB_DIRNAME)
	rm -rf lib/
	mv $(SCRATCH) lib
	rm -f $(_GUIB_DIR)/$(_GUIB_TGZNAME)

_repack_pwgui:
	$(REPACK) pwgui $(PWGUI_TARNAME) $(PWGUI_DIRNAME)
	rm pwgui.tar

_create_starkit_vfs: _clean_scratch
	cd $(PWGUI_VFS); touch make.inc make.versions; make distclean
	TOPDIR=$(HERE) PWGUI_VFS=$(PWGUI_VFS) \
	PWGUI_VERSION=$(PWGUI_VERSION) PWGUI_DIRNAME=$(PWGUI_DIRNAME) \
	GUIB_VERSION=$(_GUIB_VERSION) GUIB_DIRNAME=$(_GUIB_DIRNAME) $(STARKIT)

_message-source: 
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui-source package"
	@echo "#"
	@echo "# ------------------------------------------------------------"

_message-starkit_vfs: 
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

_clean_scratch:
	if test -d $(SCRATCH) ; then rm -rf $(SCRATCH); fi
