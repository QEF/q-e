include ../make.sys

all: libminpack pw4gwwa heada gwwa

libminpack:
	cd minpack ; $(MAKE) all

pw4gwwa: libminpack
	cd pw4gww ; $(MAKE) all

heada: pw4gwwa
	cd head ; $(MAKE) all

gwwa: libminpack
	cd gww ; $(MAKE) all

clean: 
	cd pw4gww ; $(MAKE) clean
	cd head ; $(MAKE) clean
	cd gww ; $(MAKE) clean
	cd minpack ; $(MAKE) clean
