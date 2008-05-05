
proc ::helpdoc::syntaxAppend {txt} {
    variable syntax
    
    if { [info exists syntax(count)] } {
	append syntax(txt) " "	
    } else {
	set syntax(count) 0
    }
    append syntax(txt) $txt
    incr syntax(count)    
}


proc ::helpdoc::syntaxFlush {} {
    variable syntax
    variable fid
    
    if { [info exists syntax] } {
	printf $syntax(txt)
	unset syntax
    }
}


proc ::helpdoc::manageRow {{add 0}} {
    variable rows

    set diff -1
    if { [string is integer $rows(start)] && [string is integer $rows(end)] } {
	set diff [expr $rows(end) - $rows(start)]
    }

    if { $diff > 0 } {

	# numerical arguments ...

	if { $diff < $add } {
	    return
	} elseif { $add == 2 } {
	    if { $diff > $add } { append rows(text) ". . .\n" }
	    manageRow_ $rows(end)
	} else {
	    manageRow_ [expr $rows(start) + $add]
	}
    } else {
    
	# string arguments ...
	
	if { ! [string is integer $rows(start)] } {        
	    if { $add == 0 } {
		set index $rows(start)
	    } elseif { $add < 2 } {
		set index "$rows(start)+$add"
	    } else {
		set index "$rows(end)"
		append rows(text) ". . .\n"
	    }
	    manageRow_ $index

	} elseif { ! [string is integer $rows(end)] } {

	    if { $add == 0 } {
		set index $rows(start)
	    } elseif { $add < 2 } {
		if { [string is integer $rows(start)] } {
		    set index [expr $rows(start)+$add]
		} else {
		    set index "$rows(start)+$add"
		}
	    } else {
		set index "$rows(end)"
		append rows(text) ". . .\n"
	    }
	    manageRow_ $index
	}
    }
}


proc ::helpdoc::manageRow_ {index} {
    variable rows

    foreach field $rows(line) {
	switch -- $field {
	    __conditional::begin__ {
		append rows(text) "\[ "
	    }
	    __conditional::end__ {
		append rows(text) "\] "
	    }
	    __optional::begin__ - __optional::end__ {		
		append rows(text) "$field "
	    }       
	    default {
		append rows(text) "${field}(${index}) "
	    }
	}
    }
    append rows(text) "\n"
}


proc ::helpdoc::printRows {} {
    variable rows 

    # scan $rows(text) for width

    foreach line [split $rows(text) \n] {
	
	set count 0
	
	foreach field $line {
	    
	    if { $field == "__optional::begin__" } {
		set field \{
	    }
	    if { $field == "__optional::end__" } {		
		set field \}
	    }	      

	    set len [string length $field]
	    
	    if { ! [info exists max($count)] } {
		set max($count) $len
	    } else { 
		if { $len > $max($count) } {
		    set max($count) $len
		}
	    }

	    incr count
	}
    }

    # now print

    foreach line [split $rows(text) \n] {

	set pl ""
	set count 0

	foreach field $line {
	    if { $field == "__optional::begin__" } {
		set field \{
	    }
	    if { $field == "__optional::end__" } {		
		set field \}
	    }	      

	    if { $field == "." } {
		append pl ". "
	    } else {
		append pl [format "%-$max($count)s  " $field]
	    }
	    incr count
	}

	printf ${pl}
    }
}


proc ::helpdoc::manageCol {{add 0}} {
    variable cols

    set diff -1
    if { [string is integer $cols(start)] && [string is integer $cols(end)] } {
	set diff [expr $cols(end) - $cols(start)]
    }

    if { $diff > 0 } {

	# numerical arguments ...

	if { $diff < $add } {
	    return
	} elseif { $add < 2 } {
	    lappend cols(indices) [expr $cols(start) + $add]
	} elseif { $add == 2 } {
	    if { $diff > $add } { lappend cols(indices) "..." }
	    lappend cols(indices) $cols(end)
	}
    } else {
	
	# string arguments ...
	
	if { ! [string is integer $cols(start)] } {        
	    if { $add == 0 } {
		lappend cols(indices) $cols(start)
	    } elseif { $add < 2 } {
		lappend cols(indices) "$cols(start)+$add"
	    } elseif { $add == 2 } {
		lappend cols(indices) ... 
		lappend cols(indices) $cols(end)
	    }
	} elseif { ! [string is integer $cols(end)] } {
	    
	    if { $add == 0 } {
		lappend cols(indices) $cols(start)
	    } elseif { $add < 2 } {
		lappend cols(indices) [expr $cols(start)+$add]
	    } elseif { $add == 2 } {
		lappend cols(indices) ...
		lappend cols(indices) $cols(end)
	    }
	}
    }
}


proc ::helpdoc::printCols {} {
    variable cols 

    # scan for field-width

    set extra 0

    foreach row $cols(vline) {

	switch -- $row {
	    __conditional::begin__ - __optional::begin__ {		
		incr extra 
		continue
	    }
	     __conditional::end__ - __optional::end__ {
		 continue
	     }
	}
	    
	set count 0
	
	foreach ind $cols(indices) {
	    	    
	    if { ! [info exists max($count)] } {
		set max($count) [string length ${row}(${cols(start)})]
	    }
	    
	    set _len [string length ${row}(${ind})]
	    if { $_len > $max($count) } {
		set max($count) $_len
	    }
	    incr count	    
	}	
    }

    # now print

    set ct ""
    set fie 0
    set newline 0
    foreach row $cols(vline) {
	
	if { $extra } {	    
	    switch -- $row {
		__conditional::begin__ {
		    set cbe 1
		    continue
		}
		__conditional::end__ {
		    set cen 1
		    append ct "\] "
		    continue
		}
		__optional::begin__ {		
		    set obe 1
		    continue
		}
		__optional::end__ {		
		    set oen 1
		    append ct "\} "
		    continue
		}

		default {

		    if { [info exists obe] } { incr fie } 
		    if { [info exists cbe] } { incr fie }					    
		}
	    }
	}
	 	

	if { $newline } {
	    append ct \n
	}
	append ct [::textutil::blank [expr ($extra - $fie) * 2]]
	if { [info exists obe] } { append ct "\{ " } 
	if { [info exists cbe] } { append ct "\[ " }			
	
	set count 0

	foreach ind $cols(indices) {	    
	    if { $ind == "..." } {
		append ct ". . .  "
	    } else {
		append ct [format "%-$max($count)s  " ${row}(${ind})]	    
	    }
	    incr count
	}
	
	foreach var {obe oen cbe cen} {
	    if { [info exists $var] } {
		unset $var
	    }
	}

	set fie 0
	set newline 1
    }

    # must be here, if "en" is the last row ...
    if { [info exists oen] || [info exists cen] } {
	append ct \n
    }

    printf $ct
}
