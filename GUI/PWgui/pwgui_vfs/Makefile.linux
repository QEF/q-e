# ------------------------------------------------------------------------
#
# EDIT to suit your needs ...
#
# ------------------------------------------------------------------------

PWGUI_SOURCE_TGZ =
PWGUI_VERSION = 
GUIB_VERSION  = 

# see the main.tcl.sh script for meaning of ITCL_EXACT
TK_VERSION = 8.4
ITCL_EXACT = -exact

ITCL_VERSION     = 3.2
ITK_VERSION      = 3.2
IWIDGETS_VERSION = 4.0.1

# define the library-directories to be copied
prefix=/usr/lib
ITCL_LIBDIR     = $(prefix)/itcl$(ITCL_VERSION)
ITK_LIBDIR      = $(prefix)/itk$(ITK_VERSION)
IWIDGETS_LIBDIR = $(prefix)/iwidgets$(IWIDGETS_VERSION)

# define these two vars if libitcl.so and libitk.so needs to be copied
# separately
LIBITCL_SO = /usr/lib/libitcl$(ITCL_VERSION).so.1
LIBITK_SO  = /usr/lib/libitk$(ITCL_VERSION).so.1

# if we have some dangling links or some files that are not really
# needed we should delete them ...
xdir = pwgui.vfs/lib
DELETE_LINKS = \
	$(xdir)/iwidgets$(IWIDGETS_VERSION)/demos

SDX    = $(HOME)/local/bin/sdx
TCLKIT = $(HOME)/local/bin/tclkit

# -Don't-touch-below------------------------------------------------------

PWGUI_NAME=PWgui-$(PWGUI_VERSION)
LIBDIR=pwgui.vfs/lib
HERE = $(PWD)

default: check clean tclkit-runtime main.tcl copylib starkit

from_tarball: check clean tclkit-runtime main.tcl app-pwgui copylib starkit

check:
	@if test -z "$(PWGUI_SOURCE_TGZ)" || -z "$(PWGUI_VERSION)" || -z "$(GUIB_VERSION)"; then echo; echo "*** Edit makefile first !!!"; echo; exit 1; fi

tclkit-runtime:
	if test -x $(TCLKIT) ; then cp $(TCLKIT) tclkit-runtime ; fi

main.tcl: check
	if [ ! -d pwgui.vfs ]; then MakeDir.sh $(LIBDIR); fi
	GUIB_VERSION=$(GUIB_VERSION) \
	PWGUI_VERSION=$(PWGUI_VERSION) \
	TK_VERSION=$(TK_VERSION) \
	ITCL_EXACT=$(ITCL_EXACT) \
	ITCL_VERSION=$(ITCL_VERSION) $(SHELL) main.tcl.sh

app-pwgui: check
	if test "x$(PWGUI_SOURCE_TGZ)" != "x"; then \
		cp $(PWGUI_SOURCE_TGZ) /tmp; \
		( cd /tmp; tar zxvf $(HERE)/$(PWGUI_SOURCE_TGZ) ); \
		mv /tmp/$(PWGUI_NAME) $(LIBDIR)/app-pwgui; fi

copylib:
	if test "$(ITCL_LIBDIR)" != "" ; then cp -fr $(ITCL_LIBDIR) $(LIBDIR)/; fi
	cp -fr $(ITK_LIBDIR) $(LIBDIR)/
	cp -fr $(IWIDGETS_LIBDIR) $(LIBDIR)/
	if test "$(LIBITCL_SO)" != "" ; then cp $(LIBITCL_SO) $(LIBDIR)/; fi
	if test "$(LIBITK_SO)"  != "" ; then cp $(LIBITK_SO)  $(LIBDIR)/; fi
	-rm -rf $(DELETE_LINKS)

starkit:
	$(SDX) wrap pwgui -interp $(TCLKIT) -runtime tclkit-runtime



clean:
	-rm -rf pwgui pwgui.vfs/ *~