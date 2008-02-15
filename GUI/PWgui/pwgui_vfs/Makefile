include make.sys

PWGUI_NAME=PWgui-$(PWGUI_VERSION)
VFSDIR=pwgui.vfs
LIBDIR=$(VFSDIR)/lib
HERE = $(PWD)

default: check-default clean tclkit-runtime main.tcl copylib starkit

from_tarball: check-all clean tclkit-runtime main.tcl app-pwgui copylib starkit

# check-default target: the PWGUI_SOURCE_TGZ is NOT required !!!
check-default:
	@if test -z "$(PWGUI_VERSION)" || -z "$(GUIB_VERSION)"; then echo; echo "*** Edit make.sys file first !!!"; echo; exit 1; fi

# same as check-default, but PWGUI_SOURCE_TGZ is required !!!
check-all:
	@if test -z "$(PWGUI_SOURCE_TGZ)" || -z "$(PWGUI_VERSION)" || -z "$(GUIB_VERSION)"; then echo; echo "*** Edit make.sys file first !!!"; echo; exit 1; fi

tclkit-runtime:
	if test -x $(TCLKIT) ; then cp $(TCLKIT) tclkit-runtime ; fi

main.tcl: check-default
	if test ! -d $(LIBDIR) ; then $(HERE)/MakeDir.sh $(LIBDIR); fi
	VSFDIR=$(VFSDIR) \
	GUIB_VERSION=$(GUIB_VERSION) \
	PWGUI_VERSION=$(PWGUI_VERSION) \
	TK_VERSION=$(TK_VERSION) \
	ITCL_EXACT=$(ITCL_EXACT) \
	ITCL_VERSION=$(ITCL_VERSION) $(SHELL) main.tcl.sh

app-pwgui: check-default
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
	-rm pwgui
	$(SDX) wrap pwgui -interp $(TCLKIT) -runtime tclkit-runtime



clean:
	-rm -rf pwgui *~

clean-all:
	-rm -rf pwgui $(VFSDIR) *~