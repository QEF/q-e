PW_DOC_DIR = ../../../../pwdocs

INPUT_FILES = \
	INPUT_PW \
	INPUT_PH \
	INPUT_PP \
	INPUT_CHDENS \
	INPUT_PROJWFC \
	INPUT_D3
IMAGES = \
	democritos.png \
	pwscf.png
GIF_IMAGES = $(IMAGES:.png=.gif)
LINK_FILES = \
	$(INPUT_FILES) \
	manual.tex \
	$(IMAGES)

%.pdf : %.tex
	if test ! -f $@ ; then pdflatex $<; fi
	pdflatex $<

%.gif : %.png
	convert $< $@

all: input_html compile_manual

links:
	for file in $(LINK_FILES); do ln -sf $(PW_DOC_DIR)/$$file . ; done

clean_links:
	for file in $(LINK_FILES); do \
		if test -f $$file ; then rm -f $$file ; fi ; done

input_html: links 
	./htmlise.sh PW      pw.x      INPUT_PW      > INPUT_PW.html
	./htmlise.sh PH      ph.x      INPUT_PH      > INPUT_PH.html
	./htmlise.sh PP      pw.x      INPUT_PP      > INPUT_PP.html
	./htmlise.sh ChDens  chdens.x  INPUT_CHDENS  > INPUT_CHDENS.html
	./htmlise.sh ProjWFC projwfc.x INPUT_PROJWFC > INPUT_PROJWFC.html
	./htmlise.sh D3      d3.x      INPUT_D3      > INPUT_D3.html

compile_manual: manual.pdf $(GIF_IMAGES) manual

manual.pdf: manual.tex

$(GIF_IMAGES): $(IMAGES)

manual: manual.pdf manual.tex
	rm -rf manual/	
	latex2html \
		-t "PWscf's Users Manual" \
		-html_version 3.2,math \
		-toc_depth 5 -split 5 -toc_stars -show_section_numbers \
		-local_icons \
		manual.tex
	cd manual; \
	for file in *.png; do \
		convert $$file $${file%.png}.gif; \
		rm -f $$file; \
	done; \
	for file in *.html; do \
		cp $$file /tmp/$$file; \
		cat /tmp/$$file | sed 's/\.png/\.gif/g' - | sed 's/HREF="http/NAME="http/g' - > $$file; \
		rm -f /tmp/$$file; \
	done

clean_manual:
	-rm -f manual.aux manual.log manual.out manual.toc *.gif

clean: clean_links clean_manual
	-rm -f INPUT_*.html *~
	-rm -rf manual.pdf manual/	
