# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::projwfcReadFilter {moduleObj channel} {

    # projwfc.x formatted input file should have the &INPUTPP namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel PROJWFC errMsg]
    if { $status == 0 } {	
	$moduleObj readFileWrongFormat projwfc.x $errMsg
    }
    return $channel
}