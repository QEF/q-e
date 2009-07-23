LATEX   = pdflatex
LATEX2HTML = latex2html

PDFS = constraints_HOWTO.pdf  developer_man.pdf  eps_man.pdf  user_guide.pdf
AUXS = $(PDFS:.pdf=.aux)
LOGS = $(PDFS:.pdf=.log)

all:  pdf html
pdf:  $(PDFS)
html: user_guide developer_man

$(PDFS): %.pdf: %.tex
	$(LATEX)  $<
	$(LATEX)  $<

clean:
	- rm -f $(PDFS) $(AUXS) $(LOGS) *~
	- rm -rf user_guide developer_man

developer_man: developer_man.tex
	 rm -rf developer_man/
	 latex2html \
                -t "Developer's Manual for Quantum-ESPRESSO" \
                -html_version 3.2,math \
                -toc_depth 3 -split 3 -toc_stars -show_section_numbers \
                -local_icons \
                developer_man.tex
	cd developer_man; \
           for file in *.png; do \
                convert $$file $${file%.png}.gif; \
                rm -f $$file; \
        done; \
        for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/\.png/\.gif/g' - | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' | sed 's/&nbsp;/ /g' - > $$file; \
                rm -f /tmp/$$file; \
        done

user_guide: user_guide.tex
	 rm -rf user_guide/
	 latex2html \
                -t "User's Guide for Quantum-ESPRESSO" \
                -html_version 3.2,math \
                -toc_depth 5 -split 5 -toc_stars -show_section_numbers \
                -local_icons \
                user_guide.tex
	cd user_guide; \
           for file in *.png; do \
                convert $$file $${file%.png}.gif; \
                rm -f $$file; \
        done; \
        for file in *.html; do \
                cp $$file /tmp/$$file; \
                cat /tmp/$$file | sed 's/\.png/\.gif/g' - | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' | sed 's/&nbsp;/ /g' - > $$file; \
                rm -f /tmp/$$file; \
        done

