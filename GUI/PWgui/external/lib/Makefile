
EXTERNAL_LIB_PREFIX = /home/tone/CVS/TkPWscf/external/lib
TCL_VER      = 8.3
ITCL_VER     = 3.1
IWIDGETS_VER = 3.1.0

default:
	@echo "Make: choose either links or clean"

check:

links:
	ln -s $(EXTERNAL_LIB_PREFIX)/tcl$(TCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/tk$(TCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/itcl$(ITCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/itk$(ITCL_VER) .
	ln -s $(EXTERNAL_LIB_PREFIX)/iwidgets$(IWIDGETS_VER) .

clean:
	rm -f tcl$(TCL_VER) 
	rm -f tk$(TCL_VER) 
	rm -f itcl$(ITCL_VER) 
	rm -f itk$(ITCL_VER) 
	rm -f iwidgets$(IWIDGETS_VER) 
	rm -f *~

