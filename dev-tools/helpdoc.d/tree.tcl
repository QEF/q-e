proc ::helpdoc::getFromTree {tree node key} {
    if { [$tree keyexists $node $key] } {
	return [$tree get $node $key]
    } 
    return ""
}

# TODO: reimplement as getNodeFromDescendant and getNodeFromDescendantPath
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

# TODO:
# implement ::helpdoc::getTextFromDescendantPath (this we may need
# some time). This is aka the old getDescendantText who didn't work well ...

proc ::helpdoc::getTextFromDescendant {tree node tag} {
    # PURPOSE: get text from all descendant tags named $tag
    #
    # Usage: getTextFromDescendant $tree $node $tag

    set result ""

    foreach child [$tree descendants $node] {

	set _tag [getFromTree $tree $child tag]
	
	if { $tag == $_tag } {
	    append result "[getFromTree $tree $child text] "
	}
    }    
    
    return $result
}


proc ::helpdoc::getAttributeFromDescendantPath {tree node args} {
    # PURPOSE: get the requested attribute of specified decendant
    #
    # Usage:
    #        getDescendantAttribute  $tree  $node  tag1 tag2 last_tag attribute_of_last_tag
    #
    #        where "tag1 tag2 last_tag attribute_of_last_tag" represents path to the attribute
    #

    set result ""
    set tag  [lindex $args 0]; # consider the first tag in the list of tags ...
    set att  [lindex $args end]
    
    foreach child [$tree children $node] {
	
	set _tag [getFromTree $tree $child tag]
	
	if { $tag == $_tag } {

	    # are we already at the end-path, where args = {tag attribute} ?
	    
	    if { [llength $args] == 2 } {
		# we are at the end-path, hence get the $att attribute of $_tag

		set attr [getFromTree $tree $child attributes]
		attr2array_ arr $attr
		
		if { [info exists arr($att)] } {
		    append result $arr($att)
		}
		
	    } else {
		# note yet at the end-path,
		# strip-off the current level from args and recursively re-call the proc ...
		set args1 [lrange $args 1 end]
		return [getDescendantAttribute $tree $child $args1]
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
