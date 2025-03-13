VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)

LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf dfpt_tetra.pdf dfpt_self_energy.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

HTMLS = INPUT_DYNMAT.html \
	INPUT_MATDYN.html \
	INPUT_PH.html \
	INPUT_POSTAHC.html \
	INPUT_Q2R.html

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

doc:  all
all:  defs pdf html
pdf: $(PDFS)
html: user_guide

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/
	- rm -f $(HTMLS) $(TXTS) $(XMLS) input_xx.xsl
	- rm -r $(patsubst %,../../Doc/%,$(HTMLS) $(XMLS) $(TXTS))

user_guide: user_guide.pdf
	- rm -rf user_guide/
	latex2html \
                -t "User's Guide for the PHonon package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                user_guide.tex
	cd user_guide; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** User's Guide created in user_guide/user_guide.html"
	@echo "***"
	@echo ""


defs: input_xx.xsl $(HTMLS) link_on_main_doc

$(HTMLS): %.html: %.def input_xx.xsl
	$(HELPDOC) $<

$(TXTS): %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

link_on_main_doc:
	-@( cd ../../Doc ; for file in $(HTMLS) $(XMLS) $(TXTS); do ln -fs ../PHonon/Doc/$$file . ; done )
