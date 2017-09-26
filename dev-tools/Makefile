TMP_FILES = INPUT_*.xml INPUT_*.xml.tcl *-help.tcl


pw: makedoc
	./check_gui pw

neb: makedoc
	./check_gui neb

ph: makedoc 
	./check_gui ph

pp: makedoc 
	./check_gui pp

projwfc: makedoc
	./check_gui projwfc

bands: makedoc
	./check_gui bands

dos: makedoc
	./check_gui dos

atomic: makedoc 
	./check_gui atomic

makedoc:
	(cd ..; $(MAKE) doc)

clean:
	- rm $(TMP_FILES)
