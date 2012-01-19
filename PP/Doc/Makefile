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

all: link_input_xx doc
doc: $(HTMLS) $(TXTS) $(PDFS) 

$(HTMLS): %.html: %.def
	$(HELPDOC) $<

$(TXTS): %.txt: %.def
	$(HELPDOC) $<

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

link_input_xx:
	if test -f ../../doc-def/input_xx.xsl ; then \
        (ln -s ../../doc-def/input_xx.xsl input_xx.xsl) ; \
        else \
        echo ; \
        echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
        echo ; \
        fi

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS)
	- rm -f $(HTMLS) $(TXTS) $(XMLS) input_xx.xsl *~
