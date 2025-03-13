VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)

HTMLS = INPUT_Lanczos.html \
	INPUT_Davidson.html \
	INPUT_EELS.html \
        INPUT_Magnon.html \
	INPUT_Spectrum.html

TXTS = $(HTMLS:.html=.txt)

doc:  all
all:  defs

clean:
	- rm -f input_xx.xsl
	- rm -f INPUT_*.html INPUT_*.xml INPUT_*.txt
	- rm -f $(patsubst %,../../Doc/%,$(HTMLS) $(XMLS) $(TXTS)) 

defs: input_xx.xsl $(HTMLS) $(TXTS) link_on_main_doc

input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

$(HTMLS): %.html: %.def input_xx.xsl
	$(HELPDOC) $<

$(TXTS): %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

link_on_main_doc:
	-( cd ../../Doc ; for file in $(HTMLS) $(XMLS) $(TXTS); do ln -fs ../TDDFPT/Doc/$$file . ; done )
