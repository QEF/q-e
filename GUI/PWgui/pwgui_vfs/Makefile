# ------------------------------------------------------------------------
#
# EDIT to suit your needs ...
#
# ------------------------------------------------------------------------

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

# if we have some dangling links we should delete them ...
xdir = pwgui.vfs/lib
DELETE_LINKS = \
	$(xdir)/iwidgets$(IWIDGETS_VERSION)/demos/images

SDX    = $(HOME)/local/bin/sdx
TCLKIT = $(HOME)/local/bin/tclkit

# -Don't-touch-below------------------------------------------------------

LIBDIR=pwgui.vfs/lib

all: tclkit-runtime main.tcl copylib starkit


tclkit-runtime:
	if test -x $(TCLKIT) ; then cp $(TCLKIT) tclkit-runtime ; fi

main.tcl:
	TK_VERSION=$(TK_VERSION) \
	ITCL_EXACT=$(ITCL_EXACT) \
	ITCL_VERSION=$(ITCL_VERSION) $(SHELL) main.tcl.sh

copylib:
	if test "$(ITCL_LIBDIR)" != "" ; then cp -r $(ITCL_LIBDIR) $(LIBDIR)/; fi
	cp -r $(ITK_LIBDIR) $(LIBDIR)/
	cp -r $(IWIDGETS_LIBDIR) $(LIBDIR)/
	if test "$(LIBITCL_SO)" != "" ; then cp $(LIBITCL_SO) $(LIBDIR)/; fi
	if test "$(LIBITK_SO)"  != "" ; then cp $(LIBITK_SO)  $(LIBDIR)/; fi
	-rm -f $(DELETE_LINKS)

starkit:
	$(SDX) wrap pwgui -interp $(TCLKIT) -runtime tclkit-runtime



clean:
	-rm -rf main.tcl.sh pwgui.vfs/ *~