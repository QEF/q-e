LATEX   = pdflatex

PDFS = constraints_HOWTO.pdf eps_man.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)

pdf:  $(PDFS)
all:  pdf

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	rm -f $(PDFS) $(AUXS) $(LOGS) *~