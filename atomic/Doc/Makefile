VERSION = svn
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)
LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = pseudo-gen.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

all:  pdf html defs
pdf:  $(PDFS)
html: pseudo-gen

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf pseudo-gen/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml input_xx.xsl
	- rm -rf ../../Doc/INPUT_LD1.*
	- rm -rf ../../Doc/input_xx.xsl

pseudo-gen: pseudo-gen.pdf
	 rm -rf pseudo-gen/
	 latex2html \
                -t "User's Guide for LD1" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                pseudo-gen.tex
	cd pseudo-gen; \
        for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** User's Guide created in pseudo-gen/pseudo-gen.html"
	@echo "***"
	@echo ""

defs: input_xx.xsl INPUT_LD1.html INPUT_LD1.txt link_on_main_doc

input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

INPUT_LD1.html: %.html: %.def  input_xx.xsl
	$(HELPDOC) $<
INPUT_LD1.txt: %.txt: %.def  input_xx.xsl
	$(HELPDOC) $<
link_on_main_doc:
	-@( cd ../../Doc ; ln -fs ../atomic/Doc/INPUT_LD1.html . ; \
        ln -fs ../atomic/Doc/INPUT_LD1.xml . ; \
        ln -fs ../atomic/Doc/INPUT_LD1.txt .)
