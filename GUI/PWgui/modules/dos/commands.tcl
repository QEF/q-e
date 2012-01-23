# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::dosReadFilter {moduleObj channel} {

    # dos.x formatted input file should have the &DOS namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel DOS errMsg]
    if { $status == 0 } {	
	$moduleObj readFileWrongFormat dos.x $errMsg
    }
    return $channel
}