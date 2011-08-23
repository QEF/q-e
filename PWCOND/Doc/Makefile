HELPDOC=../../dev-tools/helpdoc

doc:  all
all:  defs
clean:
	- rm -f INPUT_*.html INPUT_*.txt INPUT_*.xml
	- rm -rf ../../Doc/INPUT_PWCOND.*

defs: INPUT_PWCOND.html INPUT_PWCOND.txt link_on_main_doc
INPUT_PWCOND.html: %.html: %.def
	$(HELPDOC) $<
INPUT_PWCOND.txt: %.txt: %.def
	$(HELPDOC) $<
link_on_main_doc:
	-( cd ../../Doc ; ln -fs ../PWCOND/Doc/INPUT_PWCOND.html . ; \
	ln -fs ../PWCOND/Doc/INPUT_PWCOND.xml . ; \
	ln -fs ../PWCOND/Doc/INPUT_PWCOND.txt .)
