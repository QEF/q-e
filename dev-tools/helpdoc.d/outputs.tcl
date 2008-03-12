#
# XML
#
proc ::helpdoc::xml_escape_chr {content} {
    # replace xml special characters by escape-characters
    foreach {chr escChr} {
	& 	{\&amp;}  
	< 	{\&lt;}   
	> 	{\&gt;}   
    } {
	regsub -all -- $chr $content $escChr content
    }
    regsub -all -- '  $content {\&apos;} content
    regsub -all -- \" $content {\&quot;} content

    return $content
}


proc ::helpdoc::xml_tag_enter {tag attr content depth} {
    variable fid
    
    set indent [indent $depth]

    set sep ""
    if { $content != "" } {   
	if { [llength [split $content \n]] > 1 } {
	    set content [trimEmpty $content]
	    set sep \n
	} else {
	    set sep " "
	}
    }
    
    set content [formatString [xml_escape_chr $content]]

    if { $attr != "" } {
	puts $fid(xml) "${indent}<$tag ${attr}>${sep}${content}"
    } else {
	puts $fid(xml) "${indent}<$tag>${sep}${content}"
    }
}

proc ::helpdoc::xml_tag_leave {tag attr content depth} {
    variable fid
    puts $fid(xml) "[indent $depth]</$tag>"
}


#
# TXT
#
proc ::helpdoc::attr2array_ {arrayVar attributes} {
    upvar $arrayVar attr

    foreach {name value} [::textutil::splitx $attributes "=\"|\"\[ \n\r\\t\]|\"$"] {
	if { $name != "" } {
	    set attr($name) [string trim $value =]
	}
    }
}


proc ::helpdoc::txt_tag_enter {tag attr content depth} {
    variable txtDepth
    variable indentNum
    variable fid

    set indent  [indent $txtDepth]
    set content [formatString [trimEmpty $content]]
    
    attr2array_ arr $attr
    
    switch -exact $tag {
	namelist {
	    incr txtDepth
	    puts $fid(txt) "${indent}========================================================================\n"
	    puts $fid(txt) "${indent}[string toupper $tag] &$arr(name)"
	}
	text {
	    puts $fid(txt) [formatString $content $txtDepth]
	}
	var {
	    # if arr(type) does not exists, then var was called from vargroup
	    if { [info exists arr(type)] } {
		puts $fid(txt) ${indent}[format "%-15s %s" $arr(name) [string toupper $arr(type)]]
	    }
	}
	info {
	    puts $fid(txt) [formatString $content $txtDepth [expr 15 + 1]]
	}
	status {
	    set content "[::tclu::labelMsg {( Status} $content] )"
	    puts $fid(txt) [formatString $content $txtDepth [expr 15 + 1]]
	    #puts $fid(txt) ${indent}[format "%15s %s" {} "( Status:  $content )"]
	}
	label {
	    puts $fid(txt) "${indent}!"
	    puts $fid(txt) ${indent}[::tclu::labelMsg "${indent}! :::" [formatString $content]]
	    puts $fid(txt) "${indent}!\n" 
	}
    }
    if { $tag == "default" } {
	set content "[::tclu::labelMsg {( Default} $content] )"
	puts $fid(txt) [formatString $content $txtDepth [expr 15 + 1]]
	#puts $fid(txt) ${indent}[format "%15s %s" {} "( Default:  $content )"]
    }
}


proc ::helpdoc::txt_tag_leave {tag attr content depth} {
    variable fid 
    variable txtDepth   
    switch -exact $tag {
	namelist {
	    incr txtDepth -1
	    set indent  [indent $txtDepth]
	    puts $fid(txt) "${indent}END OF NAMELIST\n"
	}
	group - text {
	    puts $fid(txt) ""
	}
	var {
	    puts $fid(txt) "\n"
	}
    }       
}


#
# Robodoc
#
proc ::helpdoc::rbd_tag_enter {tag attr content depth} {
    variable fid 
    variable rbd_var 
    variable rbd_stack 
    variable rbd_info

    set content [formatString [trimEmpty $content]]
    
    attr2array_ arr $attr
    
    switch -exact $tag {
	input_description {
	    set rbd_stack [::struct::stack]
	    
	    set module {}
	    set rbd_info(program) unknown
	    if { [info exists arr(distribution)] } { set module $arr(distribution) }
	    if { [info exists arr(package)]      } { set module $arr(package) }
	    if { [info exists arr(program)]      } { 
		set module $module/$arr(program) 
		set rbd_info(program) $arr(program)
	    }
	    if {  $module == ""                  } { set module /input }
	    	    
	    set current_module [lindex [split $module /] end]
	    $rbd_stack push $current_module
	    
	    puts $fid(rbd) [formatString [subst {
		#****h* $module
		# DESCRIPTION
		#   Description of the input syntax for program ...
		#******
	    }]]\n
	}
	namelist {
	    set module "[$rbd_stack peek]/$arr(name)"
	    $rbd_stack push $arr(name)

	    puts $fid(rbd) [formatString [subst {
		#****n* $module
		# DESCRIPTION
		#   Description of the $arr(name) namelist.
		#******
	    }]]\n
	}
	var {
	    set name $arr(name)
	    regsub -all -- , $name + name
	    set rbd_var "#****v* [$rbd_stack peek]/$name\n"
	    append rbd_var "# NAME\n"
	    append rbd_var "#   $arr(name)\n"
	}
	info {
	    append rbd_var "# DESCRIPTION\n[::textutil::indent $content {#   }]"
	}
	status {
	    append rbd_var "# STATUS\n[::textutil::indent $content {#   }]\n"
	}
    }
    if { $tag == "default" } {
	append rbd_var "# DEFAULT\n[::textutil::indent $content {#   }]\n"
    }
}

proc ::helpdoc::rbd_tag_leave {tag attr content depth} {
    variable fid 
    variable rbd_var 
    variable rbd_stack
    switch -exact $tag {
	namelist {
	    puts $fid(rbd) "\n\# *** END of NAMELIST\n"
	    $rbd_stack pop
	}
	var {
	    puts $fid(rbd) $rbd_var
	    puts $fid(rbd) "#******\n"
	}
    }
}


proc ::helpdoc::print {tree node action} {
    variable fid
    
    set depth [$tree depth $node]

    set tag        [$tree get $node tag]
    set attributes [getFromTree $tree $node attributes]
    set content    [getFromTree $tree $node text]

    
    # XML
    
    xml_tag_${action} $tag $attributes $content $depth

    # currently disabled:

    # TXT

    #txt_tag_${action} $tag $attributes $content [expr $depth - 1]

    # robodoc

    #rbd_tag_${action} $tag $attributes $content $depth
}
