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

PDFS = constraints_HOWTO.pdf eps_man.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)

html: $(HTMLS)
pdf:  $(PDFS)
all:  html pdf

$(HTMLS): %.html: %.def
	$(HELPDOC) $<

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	rm -f $(HTMLS) $(XMLS) $(PDFS) $(AUXS) $(LOGS) *~