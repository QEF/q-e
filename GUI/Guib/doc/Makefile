# $Id: Makefile,v 1.3 2008-02-18 18:06:30 kokalj Exp $

ROBODOC      = robodoc
ROBODOC_HTML = html

%.pdf : %.tex
	if test ! -f $@ ; then pdflatex $<; fi
	pdflatex $<

default: usage.pdf usage robodoc 

robodoc:
	if test -d $(ROBODOC_HTML) ; then rm -rf $(ROBODOC_HTML) ; fi
	$(ROBODOC)

#	( cd $(ROBODOC_HTML); \
#		cat robodoc.css | sed "s/11pt/16pt/g" - > temp.css; \
#		cat temp.css    | sed "s/10pt/12pt/g" - > robodoc.css; \
#		rm -f temp.css )

usage.pdf: usage.tex

usage: usage.pdf usage.tex
	rm -rf usage/
	latex2html -split 0 usage.tex

clean:
	if test -d $(ROBODOC_HTML) ; then rm -rf $(ROBODOC_HTML); fi
	if test -d usage/ ; then rm -rf usage/; fi
	-rm -f usage.aux usage.log usage.pdf usage.toc *~
