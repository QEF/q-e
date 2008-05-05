proc ::helpdoc::tag {{level -2}} {
    # PURPOSE
    #   Return the name of the calling proc, which is used as the name
    #   of tag.
    return [namespace tail [lindex [info level $level] 0]]
}


proc helpdoc::indent {depth {extraDepth 0}} {      
    variable indentNum
    return [::textutil::blank [expr ($depth + $extraDepth) * $indentNum]]
}


proc ::helpdoc::formatString {string {depth 0}} {
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


