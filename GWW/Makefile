include ../make.sys

all: libminpack pw4gwwa heada gwwa  

libminpack:
	cd minpack ; $(MAKE) all

pw4gwwa:
	cd pw4gww ; $(MAKE) all

heada:
	cd head ; $(MAKE) all

gwwa:
	cd gww ; $(MAKE) all

clean: 
	( cd pw4gww  ; $(MAKE) $(MFLAGS) clean )
	( cd head    ; $(MAKE) $(MFLAGS) clean )
	( cd gww     ; $(MAKE) $(MFLAGS) clean )
	( cd minpack ; $(MAKE) $(MFLAGS) clean )
