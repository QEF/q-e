.SUFFIXES : .pw-in .tddfpt-in .tddfpt_pp-in .dav-in .pw-out .tddfpt-out .tddfpt_pp-out .dav-out .tddfpt-st-in .tddfpt-st-out .tddfpt-in-s2 .tddfpt-out-s2 .tddfpt-in-s2-2 .tddfpt-out-s2-2 .pw-ref .tddfpt-ref .tddfpt-ref-s2 .tddfpt-ref-s2-2 .tddfpt-eels-in .tddfpt-eels-out .pp-in .pp-out

##################################################
#   RECOMMENDED USER EDITABLE PART STARTS HERE   #
##################################################
#An example on how to use these variables are as follows
#RUNNER = openmpi-1.4.3/bin/mpirun -np 8 
#PWSCF = "../../../bin/pw.x" -ndiag 4

RUNNER =
PWSCF = "../../../bin/pw.x" 
TDDFPT_LANCZOS = "../../../bin/turbo_lanczos.x"
TDDFPT_DAVIDSON = "../../../bin/turbo_davidson.x"
TDDFPT_EELS = "../../../bin/turbo_eels.x"
TDDFPT_PP = "../../../bin/turbo_spectrum.x"
PWPP = "../../../bin/pp.x"


##################################################
#   End of recommended editable part             #
##################################################


check_pw = "../tools/check_pw.tcl"
check_tddfpt = "../tools/check_tddfpt.tcl"
check_pp = "../tools/check_pp.tcl"
check_move_cube = "../tools/check_move_cube.sh"


.pw-in.pw-out:
	$(RUNNER) $(PWSCF) < $< > $@

.tddfpt_pp-in.tddfpt_pp-out:
	$(RUNNER) $(TDDFPT_PP) < $< > $@

.tddfpt-in.tddfpt-out:
	$(RUNNER) $(TDDFPT_LANCZOS) < $< > $@

.tddfpt-st-in.tddfpt-st-out:
	$(RUNNER) $(TDDFPT_LANCZOS) < $< > $@

.tddfpt-in-s2.tddfpt-out-s2:
	$(RUNNER) $(TDDFPT_LANCZOS) < $< > $@

.tddfpt-in-s2-2.tddfpt-out-s2-2:
	$(RUNNER) $(TDDFPT_LANCZOS) < $< > $@


.dav-in.dav-out:
	$(RUNNER) $(TDDFPT_DAVIDSON) < $< > $@

.tddfpt-eels-in.tddfpt-eels-out:
	$(RUNNER) $(TDDFPT_EELS) < $< > $@

.pp-in.pp-out:
	$(RUNNER) $(PWPP) < $< > $@




