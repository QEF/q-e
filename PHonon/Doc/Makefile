VERSION = svn
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)

LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf developer_man.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)


doc:  all
all:  pdf html defs
pdf: $(PDFS)
html: user_guide developer_man

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/ developer_man/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf input_xx.xsl
	- rm -rf ../../Doc/INPUT_PH.*

user_guide: user_guide.pdf
	rm -rf user_guide/
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


developer_man: developer_man.pdf
	rm -rf developer_man/
	latex2html \
                -t "User's Guide for the PHonon package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                developer_man.tex
	cd developer_man; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** Developer manual created in developer_man/developer_man.html"
	@echo "***"
	@echo ""


defs: input_xx.xsl INPUT_PH.txt INPUT_PH.html INPUT_DYNMAT.txt INPUT_DYNMAT.html link_on_main_doc

input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

INPUT_PH.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_PH.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

INPUT_DYNMAT.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_DYNMAT.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<

link_on_main_doc:
	-( cd ../../Doc ; ln -fs ../PHonon/Doc/INPUT_PH.html . ; \
	ln -fs ../PHonon/Doc/INPUT_PH.xml . ; \
	ln -fs ../PHonon/Doc/INPUT_PH.txt . )
