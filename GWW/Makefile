include ../make.sys

#all: libminpack pw4gwwa heada gwwa
all:
	@echo "GWW no longer supported"

libminpack:
	cd minpack ; $(MAKE) all

pw4gwwa: libminpack
	cd pw4gww ; $(MAKE) all

heada: pw4gwwa
	cd head ; $(MAKE) all

gwwa: libminpack
	cd gww ; $(MAKE) all

clean: 
	cd pw4gww ; $(MAKE) $@
	cd head ; $(MAKE) $@
	cd gww ; $(MAKE) $@
	cd minpack ; $(MAKE) $@
