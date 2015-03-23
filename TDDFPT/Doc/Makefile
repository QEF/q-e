HELPDOC=../../dev-tools/helpdoc



doc:  all
all:  defs

clean:
	- rm -f input_xx.xsl
	- rm -f INPUT_*.html INPUT_*.xml
	- rm -f INPUT_{Spectrum,Lanczos,Davidson,EELS}.txt
	- rm -f ../../Doc/INPUT_{Spectrum,Lanczos,Davidson,EELS}.*


defs: link_input_xx INPUT_Lanczos.html INPUT_Lanczos.txt \
                    INPUT_Davidson.html INPUT_Davidson.txt \
                    INPUT_EELS.html INPUT_EELS.txt \
                    INPUT_Spectrum.html INPUT_Spectrum.txt \
                    link_on_main_doc

link_input_xx:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; \
	fi) ; fi)

INPUT_Lanczos.html: %.html: %.def
	$(HELPDOC) $<
INPUT_Lanczos.txt: %.txt: %.def
	$(HELPDOC) $<

INPUT_Davidson.html: %.html: %.def
	$(HELPDOC) $<
INPUT_Davidson.txt: %.txt: %.def
	$(HELPDOC) $<

INPUT_EELS.html: %.html: %.def
	$(HELPDOC) $<
INPUT_EELS.txt: %.txt: %.def
	$(HELPDOC) $<

INPUT_Spectrum.html: %.html: %.def
	$(HELPDOC) $<
INPUT_Spectrum.txt: %.txt: %.def
	$(HELPDOC) $<


link_on_main_doc:
	-@( cd ../../Doc ; \
        ln -fs ../TDDFPT/Doc/INPUT_Lanczos.html . ; \
	ln -fs ../TDDFPT/Doc/INPUT_Davidson.html . ; \
	ln -fs ../TDDFPT/Doc/INPUT_EELS.html . ; \
	ln -fs ../TDDFPT/Doc/INPUT_Spectrum.html . ; \
	ln -fs ../TDDFPT/Doc/INPUT_Lanczos.xml . ; \
	ln -fs ../TDDFPT/Doc/INPUT_Davidson.xml . ; \
	ln -fs ../TDDFPT/Doc/INPUT_EELS.xml . ; \
	ln -fs ../TDDFPT/Doc/INPUT_Spectrum.xml . ; \
        ln -fs ../TDDFPT/Doc/INPUT_Lanczos.txt . ; \
        ln -fs ../TDDFPT/Doc/INPUT_Davidson.txt . ; \
        ln -fs ../TDDFPT/Doc/INPUT_EELS.txt . ; \
	ln -fs ../TDDFPT/Doc/INPUT_Spectrum.txt .) 
