HELPDOC = ../dev-tools/helpdoc
DOC_DIR = ../Doc

HTMLS = \
	INPUT_PW.html \
	INPUT_PH.html \
	INPUT_PP.html \
	INPUT_PROJWFC.html \
	INPUT_D3.html \
	INPUT_DOS.html \
	INPUT_BANDS.html \
	INPUT_pw_export.html

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

helpdoc: $(HTMLS) $(TXTS) copy
all:  helpdoc

copy: $(HTMLS) $(TXTS)
	cp $(HTMLS) $(TXTS) $(DOC_DIR)

$(HTMLS): %.html: %.def input_xx.xsl
	$(HELPDOC) $<

$(TXTS): %.txt: %.def
	$(HELPDOC) $<

clean:
	- rm -f $(HTMLS) $(TXTS) $(XMLS) *~
