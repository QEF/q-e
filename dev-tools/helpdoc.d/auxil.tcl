proc ::helpdoc::tag {{level -2}} {
    # PURPOSE
    #   Return the name of the calling proc, which is used as the name
    #   of tag.
    return [namespace tail [lindex [info level $level] 0]]
}


proc helpdoc::indent {depth {extraDepth 0}} {
    # PURPOSE
    # return the proper number of whitespaces for the indent at level $depth
    variable indentNum
    return [::textutil::blank [expr ($depth + $extraDepth) * $indentNum]]
}


proc ::helpdoc::formatString {string {depth 0}} {
    # PURPOSE
    # return properly indented string
    variable indentNum
    set indent [indent $depth]
    return [::textutil::indent \
                [::textutil::undent \
                     [::textutil::untabify [::textutil::trimEmptyHeading $string]]] \
                $indent]
}


proc ::helpdoc::trimEmpty {text} {
    # PURPOSE 
    # Trim empty lines (this is not equal to [string trim], because the
    # beginning and ending indenation would be lost with the latter.

    regsub -- "^(\[ \t\]*\n)*" $text {} text
    regsub -- "(\[ \t\n\])*$" $text {} text
    return $text
}


proc ::helpdoc::value_of {varname} {
    # PURPOSE
    # return the value of variable or "" if the variable is not defined
    upvar $varname var
    
    if { [ info exists var] } {
	return $var
    } else {
	return ""
    }
}


proc ::helpdoc::supercardStarttag {} {
    # PURPOSE
    # return the starttag of the supercard, i.e., value of the
    # -starttag attribute if it is defined, otherwise return the supercard's ident (name)
    
    variable arr
    if { ! [array exists arr] } { return "" }

    set start [arr starttag]

    if { $start eq {} } {
	set start [arr name]
    }
    return $start    
}
