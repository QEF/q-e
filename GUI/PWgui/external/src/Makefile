########################################################################
#
# System-dependent definitions for Linux.
# Edit according to your own needs.
#
########################################################################

# what is our Tcl/Tk version ?
include ../../Make.include

#------------------------------------------------------------------------
# if you have a GNU make it is better to set the MAKE variable to point
# to gnu make
#------------------------------------------------------------------------
MAKE   = make

#------------------------------------------------------------------------
# compilers & flags
#------------------------------------------------------------------------
# static:
#--------
CFLAGS  = -static -ansi -fwritable-strings
MATH    = -lm
CC      = gcc
LDLIB   = /usr/lib/libdl.a /usr/lib/libpthread.a
# shared:
#--------
#CFLAGS  = 
#MATH    = 
#CC      = 
#LDLIB   = 


#------------------------------------------------------------------------
# X-libraries & include files
#------------------------------------------------------------------------
X_LIB     = -L/usr/X11R6/lib -lXmu -lX11 -lXext


#
# Libraries (static loading)
#
TCL_LIB   = /usr/lib/libtcl$(TCL_VERSION).a
TK_LIB    = /usr/lib/libtk$(TCL_VERSION).a
ITCL_LIB  = /usr/lib/libitcl$(ITCL_VERSION).a /usr/lib/libitclstub$(ITCL_VERSION).a
ITK_LIB   = /usr/lib/libitk$(ITCL_VERSION).a /usr/lib/libitkstub$(ITCL_VERSION).a

#
# Include directories
#
TCL_INCDIR  = -I/usr/include/tcl$(TCL_VERSION) \
	-I/usr/include/tcl$(TCL_VERSION)/tcl-private/generic
ITCL_INCDIR = -I/usr/include/tcl$(TCL_VERSION)/itcl-private/generic
TK_INCDIR   = -I/usr/include/tk$(TCL_VERSION) \
	-I/usr/include/tcl$(TCL_VERSION)/tk-private/generic
ITK_INCDIR  = -I/usr/include/tcl8.4/itk-private/generic

INCS = $(TCL_INCDIR) $(TK_INCDIR) $(ITCL_INCDIR) $(ITK_INCDIR)

.c.o:
	$(CC) $(CFLAGS) $(INCS) -c $<

default: itkwish

all: itkwish install

itkwish: tkAppInit.o
	$(CC) $(CFLAGS) tkAppInit.o -o itkwish \
		$(TK_LIB) $(TCL_LIB) $(ITK_LIB) $(ITCL_LIB) \
		$(X_LIB) $(MATH) $(LDLIB)

install:
	if test -d ../../bin; then mv itkwish ../../bin; fi

clean:
	-rm *.o 
	-rm itkwish 
	-rm *~