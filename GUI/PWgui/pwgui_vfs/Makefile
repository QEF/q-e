
include make.inc
include make.versions

VFSDIR=pwgui.vfs
LIBDIR=$(VFSDIR)/lib
HERE = $(PWD)

default: check-default clean tclkit-runtime main.tcl copylib debian_stupid starkit

check-default:
	@if test -z "$(PWGUI_VERSION)" || -z "$(GUIB_VERSION)" || -f main.tcl.sh; then echo; echo "*** Cannot use this Makefile directly; first do: \"make pwgui-starkit-vfs\" from the PWgui root directory !!!"; echo; exit 1; fi

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

copylib:
	if test "$(ITCL_LIBDIR)" != "" ; then cp -fr $(ITCL_LIBDIR) $(LIBDIR)/; fi
	cp -fr $(ITK_LIBDIR) $(LIBDIR)/
	cp -fr $(IWIDGETS_LIBDIR) $(LIBDIR)/
	if test "$(LIBITCL_SO)" != "" ; then cp $(LIBITCL_SO) $(LIBDIR)/itcl$(ITCL_VERSION)/; fi
	if test "$(LIBITK_SO)"  != "" ; then cp $(LIBITK_SO)  $(LIBDIR)/itk$(ITK_VERSION)/; fi
	if test "$(ITCL_pkg)" != "" ; then cp $(ITCL_pkg) $(LIBDIR)/itcl$(ITCL_VERSION)/; fi
	if test "$(ITK_pkg)"  != "" ; then cp $(ITK_pkg)  $(LIBDIR)/itk$(ITK_VERSION)/; fi
	-rm -rf $(DELETE_LINKS)


# in Debian the Itcl&Itk pkgIndex.tcl contain a hard-coded directory
# and the files are scattered in /usr/lib and /usr/share. Fix it.

debian_stupid:
	for file in $(LIBDIR)/itcl$(ITCL_VERSION)/pkgIndex.tcl $(LIBDIR)/itk$(ITK_VERSION)/pkgIndex.tcl; do \
	if test -f $$file; then \
		sed 's/\/usr lib x86_64-linux-gnu/$$dir/' $$file > /tmp/pkg; cp /tmp/pkg $$file; \
	fi; done 

starkit:
	-if test -f pwgui; then rm pwgui; fi
	$(SDX) wrap pwgui.kit -interp $(TCLKIT) -runtime tclkit-runtime
	$(SDX) mkpack pwgui.kit pwgui
	chmod +x pwgui
	-rm -f tclkit-runtime pwgui.kit

# -interp tclsh -runtime tclkit-runtime
#	$(SDX) mkpack oldstar.kit newstar.kit
#	$(SDX) wrap pwgui -runtime /tmp/starkit/dummy-assam
#	$(SDX) wrap pwgui -interp $(TCLKIT) -runtime tclkit-runtime

clean:
	-rm -rf pwgui pwgui.kit *~ pwgui*tgz pwgui*zip 

distclean: clean
	-rm -rf $(VFSDIR) tclkit-runtime make.inc main.tcl.sh make.versions
