HELPDOC = ../dev-tools/helpdoc
DOC_DIR = ../Doc
ATOMIC_DOC_DIR = ../atomic_doc

HTMLS = \
	INPUT_PW.html \
        INPUT_CP.html \
	INPUT_PH.html \
        INPUT_NEB.html \
	INPUT_PP.html \
	INPUT_PROJWFC.html \
	INPUT_LD1.html \
	INPUT_D3.html \
	INPUT_DOS.html \
	INPUT_BANDS.html \
	INPUT_CPPP.html \
	INPUT_PWCOND.html \
	INPUT_pw_export.html

XMLS = $(HTMLS:.html=.xml)
TXTS = $(HTMLS:.html=.txt)

helpdoc: $(HTMLS) $(TXTS) copy copy_ld1
all:  helpdoc

copy: $(HTMLS) $(TXTS)
	cp $(HTMLS) $(TXTS) $(DOC_DIR)

copy_ld1: INPUT_LD1.html INPUT_LD1.txt
	cp INPUT_LD1.html INPUT_LD1.txt $(ATOMIC_DOC_DIR)

$(HTMLS): %.html: %.def input_xx.xsl
	$(HELPDOC) $<

$(TXTS): %.txt: %.def
	$(HELPDOC) $<

clean:
	- rm -f $(HTMLS) $(TXTS) $(XMLS) *~
