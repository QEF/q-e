# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::nebReadFilter {moduleObj channel} {

    # neb.x formatted input file (neb.dat) should have the &PATH namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel PATH errMsg]
    if { $status == 0 } {	
	$moduleObj readFileWrongFormat neb.x $errMsg
    }
    return $channel
}