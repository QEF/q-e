ROBODOC      = robodoc
ROBODOC_HTML = html

%.pdf : %.tex
	if test ! -f $@ ; then pdflatex $<; fi
	pdflatex $<

default: usage.pdf usage robodoc 

robodoc:
	if test -d $(ROBODOC_HTML) ; then rm -rf $(ROBODOC_HTML) ; fi
	$(ROBODOC)

usage.pdf: usage.tex

usage: usage.pdf usage.tex
	rm -rf usage/
	latex2html -split 0 usage.tex

clean:
	if test -d $(ROBODOC_HTML) ; then rm -rf $(ROBODOC_HTML); fi
	if test -d usage/ ; then rm -rf usage/; fi
	-rm -f usage.aux usage.log usage.pdf usage.toc *~
