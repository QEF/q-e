# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::phSelectPunchFile {moduleObj variable} {

    selectFileRoot $moduleObj $variable

}


# ------------------------------------------------------------------------
# TODO: check if input file is a ph.x input file
# ------------------------------------------------------------------------
proc ::pwscf::phReadFilter {moduleObj channel} {

    # ph.x formatted input file should have the &INPUTPH namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel INPUTPH errMsg]
    if { $status == 0 } {
	$moduleObj readFileWrongFormat ph.x $errMsg
    }

    return [readFilter::default $moduleObj $channel {logical amass verbosity}]
}
