# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::phSelectPunchFile {moduleObj variable} {

    selectFileRoot $moduleObj $variable

    #variable pwscf
    #global env
    #set _dir [string trim [$moduleObj varvalue outdir] "'"]
    #
    #if { [file isdirectory $_dir] } {
    #	set dir $_dir
    #} elseif { [info exists pwscf($moduleObj,LASTDIR,punchfile)] } {
    #	set dir $pwscf($moduleObj,LASTDIR,punchfile)
    #} else {
    #	set dir $pwscf(PWD)
    #}        
    #
    #set file [tk_getOpenFile \
    #		  -initialdir $dir \
    #		  -title      "Select a Punch File"]
    #if { $file == "" } {
    #	return
    #}
    #set pwscf($moduleObj,LASTDIR,punchfile) [file dirname $file]
    #
    #set file '[file tail [file rootname $file]]'
    #$moduleObj varset $variable -value $file    
}


# ------------------------------------------------------------------------
#
# ------------------------------------------------------------------------
proc ::pwscf::phReadFilter {moduleObj channel} {

    return [readFilter::default $moduleObj $channel {logical amass}]

    ## search for the amass(X) array-element with the largest index,
    ## then varset the ntyp to the maximum found index
    #set maxIndex 0
    #set output   {}
    #while { ! [eof $channel] } {
    #	gets $channel _line	
    #	set maxIndex  [readFilter::amassIndex  $_line $maxIndex]
    #	set _line     [readFilter::logicalFlag $_line]       
    #	append output $_line\n
    #}
    #if { $maxIndex > 0 } {
    #	$moduleObj varset ntyp -value $maxIndex
    #}
    ## close the old channel
    #close $channel
    #
    ## open a new channel (i.e. temporary file) and write the 
    ## $output to it
    #set tmpfile    [::tclu::tempFile name input]
    #set newChannel [open $tmpfile w+]
    #puts  $newChannel $output
    #flush $newChannel
    #
    ## rewind the newChannel
    #seek $newChannel 0 start
    #return $newChannel
}
    
    
# # ------------------------------------------------------------------------
# #
# # ------------------------------------------------------------------------
# proc ::pwscf::phWriteFilter {moduleObj outputContent} {
#     
#     foreach line [split $outputContent \n] {
# 	if { ! [regexp {^ *ntyp} $line] } {
# 	    append output ${line}\n
# 	}
#     }
#     return $output
# }
    
    
