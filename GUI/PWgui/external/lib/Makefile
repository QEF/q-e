# what is our Tcl/Tk version ?
include ../../Make.include

default:
	@echo "Make: choose either links or clean"

check:

links:
	ln -s $(EXTERNAL_LIB_PREFIX)/tcl$(TCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/tk$(TCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/itcl$(ITCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/itk$(ITCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/iwidgets$(IWIDGETS_VER) .

copy: _copy _only_scripts

clean:
	rm -rf tcl$(TCL_VER) 
	rm -rf tk$(TCL_VER) 
	rm -rf itcl$(ITCL_VER) 
	rm -rf itk$(ITCL_VER) 
	rm -rf iwidgets$(IWIDGETS_VER) 
	rm -rf *~

_copy:
	cp -r $(EXTERNAL_LIB_PREFIX)/tcl$(TCL_VER) .
	cp -r $(EXTERNAL_LIB_PREFIX)/tk$(TCL_VER) .
	cp -r $(EXTERNAL_LIB_PREFIX)/itcl$(ITCL_VER) .
	cp -r $(EXTERNAL_LIB_PREFIX)/itk$(ITCL_VER) .
	cp -r $(EXTERNAL_LIB_PREFIX)/iwidgets$(IWIDGETS_VER) .

_only_scripts:
	if test -d tk$(TCL_VER); then \
		rm -rf tk$(TCL_VER)/demos  tk$(TCL_VER)/images; fi
	if test -d iwidgets$(IWIDGETS_VER); then \
		rm -rf iwidgets$(IWIDGETS_VER)/demos; fi