HELPDOC=../../dev-tools/helpdoc

LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)


doc:  all
all:  pdf html defs
pdf: $(PDFS)
html: user_guide

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf ../../Doc/INPUT_PW.*


user_guide: user_guide.pdf
	rm -rf user_guide/
	latex2html \
                -t "User's Guide for the PWscf package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                user_guide.tex
	cd user_guide; \
	for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** User's Guide created in user_guide/user_guide.html"
	@echo "***"
	@echo ""


defs: INPUT_PW.html INPUT_PW.txt link_on_main_doc
INPUT_PW.html: %.html: %.def
	$(HELPDOC) $<
INPUT_PW.txt: %.txt: %.def
	$(HELPDOC) $<

link_on_main_doc:
	-( cd ../../Doc ; ln -fs ../PW/Doc/INPUT_PW.html . ; \
	ln -fs ../PW/Doc/INPUT_PW.xml . ; \
	ln -fs ../PW/Doc/INPUT_PW.txt .) 
