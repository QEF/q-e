# Usage
# make        # compile all binary
# make clean  # remove ALL binaries and objects
EXAMPLE_OUT=input_pw_1 input_pw_2 input_hz_init input_hz_run out_pw_1 out_pw_2 save out_hz_init _thermo out_hz_run current_hz CRASH

all:
	@echo "Compiling scripts..."
	- cd scripts; make all; cd ..

clean : 
	cd scripts; make clean; cd ..

clean-example : 
	rm -rfv ${EXAMPLE_OUT}

veryclean: clean clean-example

move-tests:
	mkdir -pv old_test
	rm -fvr old_test/*
	mv -v ${EXAMPLE_OUT} old_test
