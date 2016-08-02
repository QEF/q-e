# makefile for pwgui
sinclude ../make.inc
# do nothing by default !!!
default:

veryclean:
	- if test -d Guib ; then ( cd Guib; $(MAKE) veryclean ) ; fi
	- if test -d PWgui ; then ( cd PWgui; $(MAKE) distclean ) ; fi
	- if test -d QE-modes ; then ( cd QE-modes; $(MAKE) veryclean ) ; fi

