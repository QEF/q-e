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
USER_GUIDES = user_guide.pdf user_guide.html autopilot_guide.pdf autopilot_guide.html
PANDOC_PDF_ENGINE_CMD = --latex-engine
PANDOC_PDF_ENGINE_CMD = --pdf-engine
PANDOC_ADDITIONAL_DEF = --metadata csquotes=yes

HTMLS = INPUT_CP.html \
	INPUT_CPPP.html

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

doc:  all
all:  pdf defs 
pdf: $(PDFS)
user_guide: $(USER_GUIDES)

user_guide.pdf: user_guide.md pandoc_template.tex tp.tex
	$(shell sed -e '/[[_TOC_]]/d' -e 's/\$$`/$$/g' -e 's/`\$$/$$/g' $< | $(PANDOC) -N --toc --metadata title="CP's user guide" --metadata customtitlepage='tp' --metadata graphics --metadata geometry='top=3cm' --metadata geometry='left=3cm' --metadata geometry='right=3cm' --metadata geometry='bottom=3cm' --metadata hyperrefoptions='linktoc=all' ${PANDOC_PDF_ENGINE_CMD}=xelatex ${PANDOC_ADDITIONAL_DEF} -o $@ --template $(word 2,$^) )

autopilot_guide.pdf: autopilot_guide.md pandoc_template.tex tp.tex
	$(shell sed -e '/[[_TOC_]]/d' -e 's/\$$`/$$/g' -e 's/`\$$/$$/g' $< | $(PANDOC) -N --toc --metadata title="CP's AUTOPILOT user guide" --metadata customtitlepage='tp' --metadata graphics --metadata geometry='top=3cm' --metadata geometry='left=3cm' --metadata geometry='right=3cm' --metadata geometry='bottom=3cm' --metadata hyperrefoptions='linktoc=all' ${PANDOC_PDF_ENGINE_CMD}=xelatex ${PANDOC_ADDITIONAL_DEF} -o $@ --template $(word 2,$^) )

user_guide.html: user_guide.md pandoc_template.html ../../Doc/quantum_espresso.png
	$(shell sed -e '/[[_TOC_]]/d' -e 's/\$$`/$$/g' -e 's/`\$$/$$/g' $< | $(PANDOC) -s --toc --metadata title="CP's user guide (v.$(VERSION))" --metadata title_image='$(word 3,$^)' --mathml -o $@ --template $(word 2,$^) )

autopilot_guide.html: autopilot_guide.md pandoc_template.html ../../Doc/quantum_espresso.png
	$(shell sed -e '/[[_TOC_]]/d' -e 's/\$$`/$$/g' -e 's/`\$$/$$/g' $< | $(PANDOC) -s --toc --metadata title="CP's AUTOPILOT user guide (v.$(VERSION))" --metadata title_image='$(word 3,$^)'  --mathml -o $@ --template $(word 2,$^) )

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) $(USER_GUIDES) *~
	- rm -rf user_guide/
	- rm -f $(HTMLS) $(TXTS) $(XMLS) input_xx.xsl
	- rm -f $(patsubst %,../../Doc/%,$(HTMLS) $(XMLS) $(TXTS))


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
	-@( cd ../../Doc ; for file in $(HTMLS) $(XMLS) $(TXTS); do ln -fs ../CPV/Doc/$$file . ; done )
