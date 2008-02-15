# ------------------------------------------------------------------------
# partial check if input file has an acceptable format
# ------------------------------------------------------------------------
proc ::pwscf::atomicReadFilter {moduleObj channel} {

    # ld1.x formatted input file should have the &INPUT namelist

    set status [::pwscf::readFilter::findNamelists $moduleObj $channel INPUT errMsg]
    if { $status == 0 } {
	$moduleObj readFileWrongFormat ld1.x $errMsg
    }
    return [readFilter::default $moduleObj $channel logical]
}
# ------------------------------------------------------------------------
# dft is a namelist variable
# We need to remove dft_ and replace dft_ with dft
# ------------------------------------------------------------------------
proc ::pwscf::atomicDFTFilter {moduleObj outputContent} { 
    set result {}

    foreach line [split $outputContent \n] {
        if { [string match {*'REPLACE_ME'*} $line] } {
            # we skip this line
        } elseif { [string match {*dft_ *} $line] } {
            # replace dft_ with dft
            
            # usage: regsub ?switches? exp string subSpec varName 
            #
            # regsub == regular-expression-substitution
            # exp     -- regular expresion to match
            # string  -- string to make the regsub 
            # subSpec -- what to do with the portion of string that matches
            #            expression (i.e. this is the replacement string)
            # varName -- where to store the result
            
            regsub -- dft_ $line dft newLine
            append result $newLine\n
        } else {
            # simply append
            append result $line\n
        }
    }
    return $result
}
#---
