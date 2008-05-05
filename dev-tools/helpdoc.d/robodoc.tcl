# Currently disabled: this file is likely to be purged in the future

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

