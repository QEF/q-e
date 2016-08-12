proc ::helpdoc::attrsToOpts_ {attrList} {
    # PURPOSE
    # Tranform attribute list to option list, i.e.:
    # {name ident type} --> {-name -ident -type}

    set optList {}
    foreach attr $attrList {
	lappend optList -$attr
    }
    return $optList
}


proc ::helpdoc::optVal2AttrVal_ {optValList} {
    # PURPOSE
    # Tranform option-value pairs to attribute value pairs, i.e.:
    # {-option1 value1 -option2 value2} --> {option1="value1" option2="value2"}

    set result ""
    foreach {opt val} $optValList {
	set attr [string trimleft $opt -]
	append result "$attr=\"$val\" "
    }
    return $result
}


proc ::helpdoc::checkIdent_ {ident} {
    # PURPOSE
    # Check if $ident is valid ident: it should not start with -, and
    # should be one word only, starting with an alphabetical
    # character"

    set ident [string trim $ident]
    set tag [tag -3]
    if { [regexp {^-} $ident] } {
	::tclu::abort "expecting ident for tag \"$tag\", but got an option $ident"
    }

    if { [llength $ident] > 1 } {
	::tclu::abort "expecting ident for tag \"$tag\" (ident should be a single word), but got a text: $ident"
    }

    if { ! [regexp {^[a-zA-Z_]} $ident] } {
	::tclu::abort "not a proper ident, $ident, for tag \"$tag\", ident start with a-z, or A-Z, or _"
    }
}

proc ::helpdoc::rootnameTag_ {args} {
    variable tree
    variable stack
    variable state
    variable elemArr
    
    set tag  [tag -2]
    set code [lindex $args end]
    set tree [::struct::tree]
    set node [$tree rootname]

    $tree set $node tag $tag   

    parseTagMsg_; puts ""

    # do tag uses ident ?
    
    #puts "tag=$tag"
    #puts "array(IDENT,*):    [array names elemArr IDENT,*]\n"
    #puts "array(ATTRLIST,*): [array names elemArr ATTRLIST,*]\n"

    if { [info exists elemArr(IDENT,$tag)] } {
	# add name="string" to attribute list
	set ident [lindex $args 0]
	checkIdent_ $ident
	set attr  "name=\"$ident\" "
	set args  [lrange $args 1 end]
    }

    # do tag use attributes ?
    
    if { [info exists elemArr(ATTRLIST,$tag)] } {
	append attr [optVal2AttrVal_ [::tclu::extractArgs \
					  [attrsToOpts_ $elemArr(ATTRLIST,$tag)]  args]]
	if { [llength $args] != 1 } {
	    # wrong attributes have been specified
	    ::tclu::abort "wrong attributes for the \"$tag\" specified, must be one of: [join $elemArr(ATTRLIST,$tag) ,]"
	}
    }

    # store attributes into the tree ...

    if { [info exists attr] } {
	$tree set $node attributes $attr    
    }

    # proceed further

    $stack push [$tree rootname]
    namespace eval tag $code
    $stack pop

    puts {[OK] - parsing finished}
}


proc ::helpdoc::elementTag_ {args} {
    variable tree
    variable stack
    variable state
    variable elemArr

    set tag  [tag -2]
    
    if { $tree == "" } {
	# an element tag has been specified before rootelement
	::tclu::abort "an element \"$tag\" specified before the rootelement \"$state(rootElem)\""
    }

    set node [$tree insert [$stack peek] end]
    set code [lindex $args end]


    $tree set $node tag $tag   

    #puts "tag=$tag"
    #puts "array(TEXT,*):     [array names elemArr TEXT,*]\n"
    #puts "array(IDENT,*):    [array names elemArr IDENT,*]\n"
    #puts "array(ATTRLIST,*): [array names elemArr ATTRLIST,*]\n"

    # do tag uses ident ?
	
    if { [info exists elemArr(IDENT,$tag)] } {
	# add name="string" to attribute list
	set name [lindex $args 0]
	parseTagMsg_ $name; 
	
	checkIdent_ $name
	set attr  "name=\"$name\" "
	set args  [lrange $args 1 end]	    
	if { $args == "" } { set code "" }
    } else {
	parseTagMsg_;
    }
    
    # do tag use attributes ?
    
    if { [info exists elemArr(ATTRLIST,$tag)] } {
	if { [llength $args] > 1 } {
	    # this is quick-and-dirty, but we need to do more cheking on order, optionality, ....
	    append attr [optVal2AttrVal_ [::tclu::extractArgs \
					      [attrsToOpts_ $elemArr(ATTRLIST,$tag)]  args]]
	    if { [llength $args] != 1 } {
		# wrong attributes have been specified
		::tclu::abort "wrong attributes for the \"$tag\" specified, must be one of: [join $elemArr(ATTRLIST,$tag) ,]"
	    }
	}
    }
	
    # TODO: checks on order, optionality, ...

    # store attributes into the tree ...
    
    if { [info exists attr] } {
	$tree set $node attributes $attr    
    }

    # we have a leaf or a complex-element ?
    
    if { [info exists elemArr(WORD,$tag)] || [info exists elemArr(STRING,$tag)] ||
	 [info exists elemArr(TEXT,$tag)] || [info exists elemArr(CLIST,$tag)] || [info exists elemArr(PLIST,$tag)] } {

	# we have a simple-element (leaf)
	$tree set $node text [lindex $args 0]
	#parseTagMsg_; puts ok
	puts ok	

    } else {
	# we have a complex element
	puts ""; # (needed for nice print-out)
	
	# proceed further

	$stack push $node
	namespace eval tag $code
	$stack pop

	parseTagMsgOK_;
    }
}


proc ::helpdoc::parseTagMsg_ {{name {}}} {
    variable tree

    set indent [uplevel 1 {indent [$tree depth $node]}]
    set tag    [string toupper [tag -3]]
    puts -nonewline "${indent}parsing $tag $name ... "    
}

proc ::helpdoc::parseTagMsgOK_ {{name {}}} {
    variable tree
    set indent [uplevel 1 {indent [$tree depth $node]}]
    set tag    [string toupper [tag -3]]
    
    if { $name == "" } {
	puts "${indent}\[OK\] - parsing $tag completed"
    } else {
	puts "${indent}\[OK\] - parsing $tag $name completed"
    }
}

