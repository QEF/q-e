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

BIN_FILES = \
	bin/itkwish

#
# per directory files ... (alphabetic order)
#
DOC_FILES = \
	doc/pwdocs/INPUT*.html \
	doc/pwdocs/users-guide.pdf \
	doc/pwdocs/users-guide/
#	doc/pwdocs/users-guide.tex 
#	doc/pwdocs/*.gif doc/pwdocs/*.png 

EXAMPLES_FILES = \
	examples/pw/*.pw.inp examples/pw/*.in \
	examples/ph/*.ph.inp \
	examples/pp/*.pp.inp \
	examples/chdens/*.ch.inp \
	examples/projwfc/*.pr.inp \
	examples/d3/*.d3.inp 

EXTERNAL_LIB_DIRS = \
	tcl$(TCL_VERSION) \
	tk$(TCL_VERSION) \
	itcl$(ITCL_VERSION) \
	itk$(ITCL_VERSION) \
	iwidgets$(IWIDGETS_VERSION) 

EXTERNAL_FILES = \
	$(addprefix external/lib/,$(EXTERNAL_LIB_DIRS)) \
	external/lib/license.terms \
	external/src/Makefile \
	external/src/tkAppInit.c \
	external/src/license.terms

IMAGES_FILES = \
	images/pwscf3.gif \
	images/pwgui-logo-small.gif \
	images/license.terms

#
# BEWARE: do not add the Guib to LIB_DIRS, as Guib should be added on the
# fly while building a package. But we need the Tcl/Tk & Itcl/Itk in
# the LIB_DIRS
#
LIB_DIRS = $(addprefix lib/,$(EXTERNAL_LIB_DIRS))

MODULES_FILES = \
	modules/*/*.tcl

SRC_FILES = \
	src/pwscf.theme \
	src/aindex \
	src/tclIndex \
	src/*.itcl

TEST_FILES = tests/run-tests.sh

PWGUI_SRC_FILES = \
	$(ROOT_FILES) $(DOC_FILES) $(EXAMPLES_FILES) \
	$(IMAGES_FILES) $(MODULES_FILES) $(SRC_FILES) $(TEST_FILES)


# nothing by default

default: 
	@cat make.usage


clean:
	-rm *~ */*~
	-rm *.tgz 
	-rm *.tar

veryclean: 
	$(MAKE) clean
	-cd doc/pwdocs; $(MAKE) clean
	-cd external/src; $(MAKE) clean
	-cd $(PWGUI_VFS); $(MAKE) clean
	-for file in `find . | egrep ~$$`; do echo "cleaning ...  $$file"; rm -f $$file; done 

distclean: veryclean
	-cd lib/; rm -f *
	-cd external/lib; $(MAKE) clean
	if test -d bin/ ; then rm -rf bin/ ; fi

_GUIB_DIR      = $(shell cd ../Guib; pwd)
_GUIB_VERSION  = $(shell GUIB=$(_GUIB_DIR); export GUIB; $(_GUIB_DIR)/guib -v)
_GUIB_TGZNAME  = guib-$(_GUIB_VERSION).tgz
_GUIB_DIRNAME  = Guib-$(_GUIB_VERSION)

cvsinit:
	-cd doc/pwdocs; $(MAKE)
	-cd lib/; rm -f Guib-*;
	-cd lib/; ln -sf $(_GUIB_DIR) $(_GUIB_DIRNAME)
	@echo ""
	@echo "NOTE: please define environmental variable PWGUI as:"
	@echo "      $(PWD)"
	@echo ""

alldist: 
	$(MAKE) pwgui-linux-i386-bin 
	$(MAKE) pwgui-source-all 
	$(MAKE) pwgui-source-notcl

pwgui-source-notcl:  \
	_message-notcl \
	_create_pwgui-notcl_tar \
	_add_guib \
	_repack_pwgui

pwgui-source-all:    \
	_message-all \
	_create_pwgui-notcl_tar \
	_add_libs \
	_add_guib \
	_repack_pwgui
	mv $(PWGUI_TARNAME).tgz $(PWGUI_TARNAME)_all.tgz

pwgui-starkit-vfs: \
	_message-starkit \
	_create_pwgui-notcl_tar \
	_add_guib \
	_create_starkit_vfs

pwgui-linux-i386-bin: \
	_message-bin \
	_create_pwgui-notcl_tar \
	_create_bin \
	_add_bin \
	_add_libs \
	_add_guib \
	_repack_pwgui
	mv $(PWGUI_TARNAME).tgz $(PWGUI_TARNAME)_linux-i386_static.tgz


_create_pwgui-notcl_tar:
	cd doc/pwdocs; $(MAKE)
	tar cvf  pwgui.tar $(PWGUI_SRC_FILES)

_add_bin:
	tar rvf  pwgui.tar $(BIN_FILES)

_create_bin:
	if test ! -d bin/ ; then mkdir bin ; fi
	cd external/src; $(MAKE) all

_add_libs:
	tar rvfh  pwgui.tar $(EXTERNAL_FILES) 
	tar rvf   pwgui.tar $(LIB_DIRS)

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
	TOPDIR=$(HERE) PWGUI_VFS=$(PWGUI_VFS) \
	PWGUI_VERSION=$(PWGUI_VERSION) PWGUI_DIRNAME=$(PWGUI_DIRNAME) \
	GUIB_VERSION=$(_GUIB_VERSION) GUIB_DIRNAME=$(_GUIB_DIRNAME) $(STARKIT)

_clean_scratch:
	if test -d $(SCRATCH) ; then rm -rf $(SCRATCH); fi

_message-notcl: _check_guib
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui-source-notcl package"
	@echo "#"
	@echo "# ------------------------------------------------------------"

_message-all: _check_guib _check_lib-links
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui-source-all package"
	@echo "#"
	@echo "# ------------------------------------------------------------"

_message-starkit: _check_guib
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui strakit (StandAlone Runtime tclKIT) VFS"
	@echo "#"
	@echo "# ------------------------------------------------------------"

_message-bin: _check_guib _check_lib-links
	@echo ""
	@echo "# ------------------------------------------------------------"
	@echo "#"
	@echo "# Creating pwgui-linux-i386-bin package"
	@echo "#"
	@echo "# ------------------------------------------------------------"

_check_guib:
	@if test -z "$(_GUIB_VERSION)" ; then \
		echo ""; echo "WARNING: GUIB package not found ..."; echo ""; exit 1; \
	fi	

_check_lib-links:
	cd lib;	\
	for link in $(EXTERNAL_LIB_DIRS); do \
		if test ! -e ../external/lib/$$link; then \
		echo ""; \
		echo "WARNING: cannot create lib-links ..."; \
		echo "         the creation of self-contained packages not possible ..."; echo ""; exit 1; \
		else ln -sf ../external/lib/$$link . ; fi \
	done

_clean_lib-links:
	cd lib; for link in $(EXTERNAL_LIB_DIRS); do rm -f $$link ; done
