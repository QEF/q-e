# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::chdensReadFilter {moduleObj channel} {

    # chdens.x formatted input file should have the &INPUT namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel INPUT errMsg]
    if { $status == 0 } {	
	$moduleObj readFileWrongFormat chdens.x $errMsg
    }
    return $channel
}