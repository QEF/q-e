########################################################################
#
# System-dependent definitions for Linux.
# Edit according to your own needs.
#
########################################################################

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
TCL_LIB   = /usr/lib/libtcl8.3.a
ITCL_LIB  = /usr/lib/libtk8.3.a
TK_LIB    = /usr/lib/libitcl3.1.a
ITK_LIB   = /usr/lib/libitk3.1.a

#
# Include directories
#
TCL_INCDIR  = -I/usr/include/tcl8.3 -I/usr/include/tcl8.3/tcl-private/generic
ITCL_INCDIR = -I/usr/include/tcl8.3/itcl-private/generic
TK_INCDIR   = -I/usr/include/tk8.3 -I/usr/include/tcl8.3/tk-private/generic
ITK_INCDIR  = -I/usr/include/tcl8.3/itk-private/generic

INCS = $(TCL_INCDIR) $(TK_INCDIR) $(ITCL_INCDIR) $(ITK_INCDIR)

.c.o:
	$(CC) $(CFLAGS) $(INCS) -c $<

all: tkAppInit.o
	$(CC) $(CFLAGS) tkAppInit.o -o itkwish \
		$(TK_LIB) $(TCL_LIB) $(ITK_LIB) $(ITCL_LIB) \
		$(X_LIB) $(MATH) $(LDLIB)

clean:
	-rm *.o itkwish *~