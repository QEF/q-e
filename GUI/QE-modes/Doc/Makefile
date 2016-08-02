LATEX      = pdflatex
LATEX2HTML = latex2html

VERSION = svn

PDFS = user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

doc: all
all: version pdf html
pdf: $(PDFS)
html: user_guide

version:
	echo '\\def\\version{$(VERSION)}' > version.tex

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

user_guide: user_guide.pdf
	- rm -rf user_guide/
	latex2html \
		-t "User's Guide for the QE-modes package" \
                -html_version 3.2,math \
                -toc_depth 5 -split 0 -toc_stars -show_section_numbers \
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

clean:
	- rm -f $(AUXS) $(LOGS) $(OUTS) $(TOCS)

veryclean: clean
	- rm $(PDFS)
	- rm version.tex
	- rm *~
	- rm -rf user_guide/

distclean: veryclean
