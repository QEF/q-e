# $Id: Makefile,v 1.1 2004-02-18 11:29:17 kokalj Exp $

ROBODOC = robodoc

%.pdf : %.tex
	if test ! -f $@ ; then pdflatex $<; fi
	pdflatex $<

default: usage.pdf usage robodoc 

robodoc:
	$(ROBODOC)
	( mkdir html; cd html; \
		cat robodoc.css | sed "s/11pt/16pt/g" - > temp.css; \
		cat temp.css    | sed "s/10pt/12pt/g" - > robodoc.css; \
		rm -f temp.css )

usage.pdf: usage.tex

usage: usage.pdf usage.tex
	rm -rf usage/
	latex2html -split 0 usage.tex

clean:
	if test -d html/  ; then rm -rf html/; fi
	if test -d usage/ ; then rm -rf usage/; fi
	-rm -f usage.aux usage.log usage.pdf usage.toc *~
