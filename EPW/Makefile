# Makefile for EPW
sinclude ../make.inc

# see comments in file tldeps for explanations
TLDEPS=phlibs libw90 pp
include ../install/tldeps

PSEUDO_DIR = ../pseudo

CD_PSEUDO = $(PSEUDO_DIR)/Cd_ONCV_PBE_SR.upf
FE_PSEUDO = $(PSEUDO_DIR)/Fe.pz.upf
O_PSEUDO = $(PSEUDO_DIR)/O_ONCV_PBE_SR.upf

CD_PSEUDO_URL = https://gitlab.com/epw/q-e/-/raw/EPW-6.0s/pseudo/Cd_ONCV_PBE_SR.upf
FE_PSEUDO_URL = https://gitlab.com/epw/q-e/-/raw/EPW-6.0s/pseudo/Fe.pz.upf
O_PSEUDO_URL = https://gitlab.com/epw/q-e/-/raw/EPW-6.0s/pseudo/O_ONCV_PBE_SR.upf

.PHONY: all epw clean distclean release download_pseudos

all: epw

epw: download_pseudos
	(cd src ; make all )
	- (mkdir -p bin)
	(cd bin ; ln -fs ../src/epw.x . )
	(cd ZG/src ; make )

download_pseudos: $(CD_PSEUDO) $(FE_PSEUDO) $(O_PSEUDO)

$(PSEUDO_DIR):
	mkdir -p $(PSEUDO_DIR)

$(CD_PSEUDO): | $(PSEUDO_DIR)
	@echo "Downloading Cd_ONCV_PBE_SR.upf..."
	@if command -v curl >/dev/null 2>&1; then \
		curl -L -o $@ $(CD_PSEUDO_URL); \
	elif command -v wget >/dev/null 2>&1; then \
		wget -O $@ $(CD_PSEUDO_URL); \
	else \
		echo "Error: neither curl nor wget is available."; \
		exit 1; \
	fi

$(FE_PSEUDO): | $(PSEUDO_DIR)
	@echo "Downloading Fe.pz.upf..."
	@if command -v curl >/dev/null 2>&1; then \
		curl -L -o $@ $(FE_PSEUDO_URL); \
	elif command -v wget >/dev/null 2>&1; then \
		wget -O $@ $(FE_PSEUDO_URL); \
	else \
		echo "Error: neither curl nor wget is available."; \
		exit 1; \
	fi

$(O_PSEUDO): | $(PSEUDO_DIR)
	@echo "Downloading O_ONCV_PBE_SR.upf..."
	@if command -v curl >/dev/null 2>&1; then \
        curl -L -o $@ $(O_PSEUDO_URL); \
    elif command -v wget >/dev/null 2>&1; then \
        wget -O $@ $(O_PSEUDO_URL); \
    else \
        echo "Error: neither curl nor wget is available."; \
        exit 1; \
    fi

clean:
	(cd src ; make clean )
	(cd ZG/src ; make clean )

distclean : clean
	rm -f src/epw.x bin/epw.x
	rm -f src/nscf2supercond.x bin/nscf2supercond.x
	rm -f ZG/src/*.x
	rm -f $(CD_PSEUDO) $(FE_PSEUDO) $(O_PSEUDO)

release:
	cd ../ ; cp -r EPW EPW-release; cd EPW-release ; \
	rm -f src/*.o src/*.mod src/*.F90 src/*~ ; \
	rm bin/*.x ; \
	rm -rf examples/*/epw/out/* examples/*/epw/tmp/* \
	 examples/*/phonons/out/* examples/*/phonons/tmp/* \
	 examples/*/phonons/save/* ; \
	cd .. ; tar cfz EPW/EPW-release.tgz EPW-release ; \
	rm -rf EPW-release ; cd EPW
