VERSION := $(shell awk -F\' '/version_number/{print $$2}' ../../include/qe_version.h )
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)
LATEX   = pdflatex
LATEX2HTML = latex2html
PANDOC = pandoc

PDFS = 
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)
USER_GUIDES = user_guide.pdf user_guide.html
PANDOC_PDF_ENGINE_CMD = --latex-engine
PANDOC_PDF_ENGINE_CMD = --pdf-engine
PANDOC_ADDITIONAL_DEF = --metadata csquotes=yes

doc:  all
all:  pdf defs 
pdf: $(PDFS)
user_guide: $(USER_GUIDES)

user_guide.pdf: user_guide.md pandoc_template.tex tp.tex
	$(shell sed -e '/[[_TOC_]]/d' -e 's/\$$`/$$/g' -e 's/`\$$/$$/g' $< | $(PANDOC) -N --toc --metadata title='CP user guide' --metadata customtitlepage='tp' --metadata graphics --metadata geometry='top=3cm' --metadata geometry='left=3cm' --metadata geometry='right=3cm' --metadata geometry='bottom=3cm' --metadata hyperrefoptions='linktoc=all' ${PANDOC_PDF_ENGINE_CMD}=xelatex ${PANDOC_ADDITIONAL_DEF} -o $@ --template $(word 2,$^) )

user_guide.html: user_guide.md
	$(shell sed -e '/[[_TOC_]]/d' -e 's/\$$`/$$/g' -e 's/`\$$/$$/g' $< | $(PANDOC) -s --toc --metadata title='CP user guide' --mathml -o $@ )

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) $(USER_GUIDES) *~
	- rm -rf user_guide/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml input_xx.xsl
	- rm -rf ../../Doc/INPUT_CP*.*
	- rm -rf ../../Doc/input_xx.xsl


defs: input_xx.xsl INPUT_CP.html INPUT_CP.txt INPUT_CPPP.html INPUT_CPPP.txt link_on_main_doc
INPUT_CP.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_CP.txt: %.txt: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_CPPP.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_CPPP.txt: %.txt: %.def input_xx.xsl
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
	-@( cd ../../Doc ; ln -fs ../CPV/Doc/INPUT_CP.html . ; \
	ln -fs ../CPV/Doc/INPUT_CP.xml . ; \
	ln -fs ../CPV/Doc/INPUT_CP.txt . ; \
	ln -fs ../CPV/Doc/INPUT_CPPP.html . ; \
	ln -fs ../CPV/Doc/INPUT_CPPP.xml . ; \
	ln -fs ../CPV/Doc/INPUT_CPPP.txt)
