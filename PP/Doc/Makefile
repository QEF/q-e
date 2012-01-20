HELPDOC = ../../dev-tools/helpdoc
LATEX   = pdflatex

HTMLS = INPUT_PP.html \
	INPUT_PROJWFC.html \
	INPUT_DOS.html \
	INPUT_BANDS.html \
	INPUT_pw_export.html

PDFS = eps_man.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

doc: all
all: pdf defs 
pdf: $(PDFS)

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

defs: link_input_xx $(HTMLS) link_on_main_doc

$(HTMLS): %.html: %.def
	$(HELPDOC) $<

$(TXTS): %.txt: %.def
	$(HELPDOC) $<

link_input_xx:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
        (ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
        else \
        echo ; \
        echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
        echo ; \
        fi) ; fi)

link_on_main_doc:
	-@( cd ../../Doc ; ln -fs ../PP/Doc/INPUT_PP.html . ; \
	ln -fs ../PP/Doc/INPUT_PP.xml . ; \
	ln -fs ../PP/Doc/INPUT_PP.txt . ; \
	ln -fs ../PP/Doc/INPUT_PROJWFC.html . ; \
	ln -fs ../PP/Doc/INPUT_PROJWFC.xml . ; \
	ln -fs ../PP/Doc/INPUT_PROJWFC.txt . ; \
	ln -fs ../PP/Doc/INPUT_DOS.html . ; \
	ln -fs ../PP/Doc/INPUT_DOS.xml . ; \
	ln -fs ../PP/Doc/INPUT_DOS.txt . ; \
	ln -fs ../PP/Doc/INPUT_BANDS.html . ; \
	ln -fs ../PP/Doc/INPUT_BANDS.xml . ; \
	ln -fs ../PP/Doc/INPUT_BANDS.txt . ; \
	ln -fs ../PP/Doc/INPUT_pw_export.html . ; \
	ln -fs ../PP/Doc/INPUT_pw_export.xml . ; \
	ln -fs ../PP/Doc/INPUT_pw_export.txt .)

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS)
	- rm -f $(HTMLS) $(TXTS) $(XMLS) input_xx.xsl *~
	- rm -rf ../../Doc/$(HTMLS)
	- rm -rf ../../Doc/$(XMLS)
	- rm -rf ../../Doc/$(TXTS)
