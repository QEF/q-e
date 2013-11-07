LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = constraints_HOWTO.pdf  developer_man.pdf  user_guide.pdf plumed_quick_ref.pdf brillouin_zones.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)
OUTS = $(PDFS:.pdf=.out)
TOCS = $(PDFS:.pdf=.toc)

all:  pdf html
pdf:  $(PDFS)
html: user_guide developer_man

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) $(OUTS) $(TOCS) *~
	- rm -rf user_guide/ developer_man/
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml qe-input-ref.html

developer_man: developer_man.pdf
	rm -rf developer_man/
	latex2html \
                -t "Developer's Manual for Quantum-ESPRESSO" \
                -html_version 3.2,math \
                -toc_depth 3 -split 3 -toc_stars -show_section_numbers \
                -local_icons -image_type png \
                developer_man.tex
	cd developer_man; \
        for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' - > $$file; \
                rm -f /tmp/$$file; \
        done
	@echo ""
	@echo "***"
	@echo "*** Developer's Manual created in developers_manual/developers_manual.html"
	@echo "***"
	@echo ""

user_guide: user_guide.pdf
	 rm -rf user_guide/
	 latex2html \
                -t "User's Guide for Quantum-ESPRESSO" \
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

brillouin_zones: brillouin_zones.pdf
