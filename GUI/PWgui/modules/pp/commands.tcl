proc ::pwscf::ppReadFilter {moduleObj channel} {

    # pp.x formatted input file should have the &INPUTPP namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel INPUTPP errMsg]
    if { $status == 0 } {
	$moduleObj readFileWrongFormat pp.x $errMsg
    }
    return [readFilter::default $moduleObj $channel logical]
}
