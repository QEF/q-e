#
# TXT
#

proc ::helpdoc::attr2array_ {arrayName attributes} {
    upvar $arrayName arrayVar

    catch {array unset arrayVar}; # EXPERIMENTAL: this should be the
			          # desired behavior, because one wants to
			          # tranform attribute-list to associative
			          # array, hence the previous key-value
			          # pairs should be cleared
    
    foreach {name value} [::textutil::splitx $attributes "=\"|\"\[ \n\r\\t\]|\"$"] {
	if { $name != "" } {
	    set arrayVar($name) [string trim $value =]
	}
    }
}

proc ::helpdoc::arr {elem} {
    variable arr

    if { [info exists arr($elem)] } {
	return $arr($elem)
    } 
    return ""
}

proc ::helpdoc::printf {content {extraSpace 0}} {
    variable txtDepth
    variable indentNum
    variable fid

    set indent [indent $txtDepth]
    if { $extraSpace > 0 } {
	set indent $indent[::textutil::blank $extraSpace]
    }
    foreach line [split $content \n] {
	puts $fid(txt) ${indent}$line
    }
}

proc helpdoc::printfNormalize {content} {
    variable txtDepth
    variable fid
    
    puts $fid(txt) [formatString $content $txtDepth]
}


proc helpdoc::labelMsg {label msg} {
    set il 1
    set len [string length $label]
    set message {}
    foreach line [split [string trim $msg] \n] {
        if { $il == 1 } {
            append message [::format "%${len}s %s" $label $line]
            incr il
        } else {
            append message [::format "\n%${len}s %s" {} $line]
        }
    }
    return $message
}

proc ::helpdoc::txt_ref_link {content} {
    set re_ref  {(@ref)\s+(\w+([%]\w)*)}
    set re_link {(@link)\s+([.,;:]*[\w\+-]+([.,;:][\w\+-]+)*)}
    set re "($re_ref|$re_link)"
    return [regsub -all $re $content {"\3"}]
}
proc ::helpdoc::txt_tag_enter {tree node tag attr content depth} {
    variable txtDepth
    variable indentNum
    variable fid
    variable arr
    variable vargroup
    variable dimensiongroup
    variable colgroup
    variable rowgroup
    variable card
    variable mode
    variable rows
    variable cols
    variable info
    variable options
    variable options_first
    
    if { [info exists arr] } {
	unset arr
    }

    set content [formatString [trimEmpty [txt_atTags [txt_ref_link $content]]]]
    #set content [formatString [trimEmpty  $content]]
    attr2array_ arr $attr

    global sourcedir
    source [file join $sourcedir txt_enter.tcl]
}


proc ::helpdoc::txt_tag_leave {tree node tag attr content depth} {
    variable fid 
    variable txtDepth   
    variable vargroup
    variable dimensiongroup
    variable colgroup
    variable rowgroup
    variable mode
    variable card
    variable rows
    variable cols
    variable arr
    variable options
    variable options_first

    attr2array_ arr $attr
    global sourcedir
    source [file join $sourcedir txt_leave.tcl]
}


proc ::helpdoc::txt_subtree {tree node newMode} {
    variable mode

    lappend mode $newMode

    set newTree [::struct::tree]
    $newTree deserialize [$tree serialize $node]

    $newTree walkproc [$newTree rootname] -order both txt_subtree_print
    $newTree destroy	

    ::tclu::lpop mode
}


proc ::helpdoc::txt_subtree_print {tree node action} {
    set depth [$tree depth $node]

    set tag        [$tree get $node tag]
    set attributes [getFromTree $tree $node attributes]
    set content    [getFromTree $tree $node text]
    
    txt_tag_${action} $tree $node $tag $attributes $content [expr $depth - 1]
}


proc ::helpdoc::printableVarDescription {tree node} {
    variable mode

    # Purpose: the description of variable in the card is printed only
    # when at least one of info, status or see records is present.

    set Info   [getTextFromDescendant $tree $node info]
    set Status [getTextFromDescendant $tree $node status]
    set See    [getTextFromDescendant $tree $node see]
    set Opt    [getTextFromDescendant $tree $node opt]

    if { ! [::tclu::lpresent $mode card] || ($Info != "" || $Status != "" || $See != "" || $Opt != "") } {
	return 1
    } 

    return 0
}
