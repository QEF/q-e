# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::hpReadFilter {moduleObj channel} {

    # hp.x formatted input file should have the &INPUTHP namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel INPUTHP errMsg]
    if { $status == 0 } {	
	$moduleObj readFileWrongFormat hp.x $errMsg
    }
    return [readFilter::default $moduleObj $channel {logical {ntyp {skip_type equiv_type perturb_only_atom}}}]
}
