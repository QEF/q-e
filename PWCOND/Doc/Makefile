VERSION = svn
HELPDOC = ../../dev-tools/helpdoc -version $(VERSION)

doc:  all
all:  defs
clean:
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf input_xx.xsl
	- rm -rf ../../Doc/INPUT_PWCOND.*

defs: input_xx.xsl INPUT_PWCOND.html INPUT_PWCOND.txt link_on_main_doc
input_xx.xsl:
	@(if test ! -f input_xx.xsl; then \
	(if test -f ../../dev-tools/input_xx.xsl; then \
	(ln -sf ../../dev-tools/input_xx.xsl input_xx.xsl) ; \
	else \
	echo ; \
	echo "  Sorry, can not find input_xx.xsl html style file !!!" ; \
	echo ; exit 1 ; \
	fi) ; fi)

INPUT_PWCOND.html: %.html: %.def input_xx.xsl
	$(HELPDOC) $<
INPUT_PWCOND.txt: %.txt: %.def  input_xx.xsl
	$(HELPDOC) $<
link_on_main_doc:
	-@( cd ../../Doc ; ln -fs ../PWCOND/Doc/INPUT_PWCOND.html . ; \
	ln -fs ../PWCOND/Doc/INPUT_PWCOND.xml . ; \
	ln -fs ../PWCOND/Doc/INPUT_PWCOND.txt .)
