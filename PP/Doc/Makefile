VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)
LATEX   = pdflatex
LATEX2HTML = latex2html


HTMLS = INPUT_BAND_INTERPOLATION.html \
	INPUT_BANDS.html \
	INPUT_bgw2pw.html \
	INPUT_D3HESS.html \
	INPUT_DOS.html \
	INPUT_IMPORTEXPORT_BINARY.html \
	INPUT_molecularpdos.html \
	INPUT_OSCDFT_ET.html \
	INPUT_OSCDFT_PP.html \
	INPUT_PPACF.html \
	INPUT_PP.html \
	INPUT_PPRISM.html \
	INPUT_PROJWFC.html \
	INPUT_pw2bgw.html \
	INPUT_pw2gw.html \
	INPUT_pw2wannier90.html

PDFS = eps_man.pdf user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

doc: all
all: defs pdf html
pdf: $(PDFS)
html: user_guide

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

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
	-@( cd ../../Doc ; for file in $(HTMLS) $(XMLS) $(TXTS); do ln -fs ../PP/Doc/$$file . ; done )

user_guide: user_guide.pdf
	- rm -rf user_guide/
	latex2html \
		-t "User's Guide for the PP package" \
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

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS)
	- rm -f $(HTMLS) $(TXTS) $(XMLS) input_xx.xsl *~
	- rm -f $(patsubst %,../../Doc/%,$(HTMLS) $(XMLS) $(TXTS))
