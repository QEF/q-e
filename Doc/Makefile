LATEX   = pdflatex

PDFS = constraints_HOWTO.pdf  developer_man.pdf  eps_man.pdf  user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)

pdf:  $(PDFS)
all:  pdf

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) *~