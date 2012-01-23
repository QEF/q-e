# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::bandsReadFilter {moduleObj channel} {

    # bands.x formatted input file should have the &BANDS namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel BANDS errMsg]
    if { $status == 0 } {	
	$moduleObj readFileWrongFormat bands.x $errMsg
    }
    return $channel
}