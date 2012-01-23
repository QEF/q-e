TMP_FILES = INPUT_*.xml INPUT_*.xml.tcl *-help.tcl


pw: makedoc
	./check_gui pw

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

d3: makedoc 
	./check_gui d3

makedoc:
	(cd ..; $(MAKE) doc)

clean:
	- rm $(TMP_FILES)