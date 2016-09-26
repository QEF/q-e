VERSION = svn
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)
LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

doc:  all
all:  pdf html defs
pdf:  $(PDFS)
html: user_guide

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf input_xx.xsl
	- rm -rf ../../Doc/INPUT_NEB.*

user_guide: user_guide.pdf
	 rm -rf user_guide/
	 latex2html \
                -t "User's Guide for The Quantum ESPRESSO Nudged Elastic Band" \
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

defs: input_xx.xsl INPUT_NEB.html INPUT_NEB.txt link_on_main_doc
input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

INPUT_NEB.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_NEB.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<
link_on_main_doc:
	-@( cd ../../Doc ; ln -fs ../NEB/Doc/INPUT_NEB.html . ; \
	ln -fs ../NEB/Doc/INPUT_NEB.xml . ; \
	ln -fs ../NEB/Doc/INPUT_NEB.txt .)
