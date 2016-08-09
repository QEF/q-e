proc ::helpdoc::getFromTree {tree node key} {
    if { [$tree keyexists $node $key] } {
	return [$tree get $node $key]
    } 
    return ""
}


proc ::helpdoc::getDescendantNodes {tree node args} {
    # Usage: getDescendantNodes $tree $node tag1 tag2 last_tag
    # get all descendant node's pointers that matches

    set result ""
    set tag [lindex $args 0]

    foreach child [$tree children $node] {

	set _tag [getFromTree $tree $child tag]
	
	if { $tag == $_tag } {
	    if { $tag == $args } {
		append result "$child "
	    } else {
		set args1 [lrange $args 1 end]
		return [getDescendantNodes $tree $child $args1]
	    }
	}
    }
    
    return $result
}


proc ::helpdoc::getDescendantText {tree node args} {
    # Usage: getDescendantText $tree $node tag1 tag2 last_tag
    # Beware: it will get the text from all tags that matches

    set result ""
    set tag [lindex $args 0]

    foreach child [$tree children $node] {

	set _tag [getFromTree $tree $child tag]
	
	if { $tag == $_tag } {
	    if { $tag == $args } {
		append result "[getFromTree $tree $child text] "
	    } else {
		set args1 [lrange $args 1 end]
		return [getDescendantText $tree $child $args1]
	    }
	}
    }
    
    return $result
}


proc ::helpdoc::getDescendantAttribute {tree node args} {

    # Usage: getDescendantText $tree $node tag1 tag2 last_tag attribute_of_last_tag
    # Beware: it will get the requested attribute from all tags that matches

    set result ""
    set tag [lindex $args 0]
    set att [lindex $args end]
    
    foreach child [$tree children $node] {
	
	set _tag [getFromTree $tree $child tag]
	
	if { $tag == $_tag } {
	    
	    if { [llength $args] == 2 } {
		
		# ok _tag is the attribute

		set attr [getFromTree $tree $child attributes]
		attr2array_ arr $attr
		
		if { [info exists arr($att)] } {
		    append result $arr($att)
		}
		
	    } else {
		set args1 [lrange $args 1 end]
		return [getDescendantAttribute $tree $child $args1]
	    }
	}
    }
    
    return $result
}
