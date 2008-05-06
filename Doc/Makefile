HELPDOC = ../dev-tools/helpdoc
LATEX   = pdflatex

HTMLS = \
	INPUT_PW.html \
	INPUT_PH.html \
	INPUT_PP.html \
	INPUT_PROJWFC.html \
	INPUT_LD1.html \
	INPUT_D3.html

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

PDFS = constraints_HOWTO.pdf eps_man.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)

helpdoc: $(HTMLS) $(TXTS)
pdf:  $(PDFS)
all:  helpdoc pdf

$(HTMLS): %.html: %.def
	$(HELPDOC) $<

$(TXTS): %.txt: %.def
	$(HELPDOC) $<

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	rm -f $(HTMLS) $(TXTS) $(XMLS) $(PDFS) $(AUXS) $(LOGS) *~