# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::d3ReadFilter {moduleObj channel} {

    # d3.x formatted input file should have the &INPUTPH namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel INPUTPH errMsg]
    if { $status == 0 } {
	$moduleObj readFileWrongFormat d3.x $errMsg
    }
    return [readFilter::default $moduleObj $channel {logical amass}]
}