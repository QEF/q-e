TMP_FILES = INPUT_*.xml INPUT_*.xml.tcl *-help.tcl


pw: makedoc
	./check_gui pw

ph: makedoc 
	./check_gui ph

pp: makedoc 
	./check_gui pp

projwfc: makedoc
	./check_gui projwfc

atomic: makedoc 
	./check_gui atomic

d3: makedoc 
	./check_gui d3

makedoc:
	(cd ../doc-def; $(MAKE) helpdoc)

clean:
	- rm $(TMP_FILES)