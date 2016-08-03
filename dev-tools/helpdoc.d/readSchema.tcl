
namespace eval ::helpdoc::schema {

    # here is the definition of Tcl-commands that are used in schema

    proc rootelement {name code} { uplevel 1 [list ::helpdoc::rootelement $name $code] }
    proc element     {name code} { uplevel 1 [list ::helpdoc::element $name $code] }
    proc attribute   {name code} { uplevel 1 [list ::helpdoc::attribute $name $code] }    
    proc define      {name code} { uplevel 1 [list ::helpdoc::define $name $code] }    
    proc text        {}          { uplevel 1 [list ::helpdoc::text] }
    proc string      {}          { uplevel 1 [list ::helpdoc::String] }
    proc clist       {}          { uplevel 1 [list ::helpdoc::clist] }
    proc ref         {name}      { uplevel 1 [list ::helpdoc::ref $name] }
    proc ident       {}          { uplevel 1 [list ::helpdoc::ident] }
    proc optional    {code}      { uplevel 1 [list ::helpdoc::optional $code] }    
    proc interleave  {code}      { uplevel 1 [list ::helpdoc::interleave $code] }
    proc choice      {code}      { uplevel 1 [list ::helpdoc::choice $code] }
    proc ancestorElements {}     { uplevel 1 [list ::helpdoc::ancestorElements] }		      
    proc ?           {code}      { uplevel 1 [list ::helpdoc::? $code] }
    proc *           {code}      { uplevel 1 [list ::helpdoc::* $code] }
    proc +           {code}      { uplevel 1 [list ::helpdoc::+ $code] } 
}

# actual implementation of commands ...

proc ::helpdoc::rootelement {name code} {
    variable elemList
    variable itemList
    variable stackArr
    variable state

    parseMsg_ $name; puts ""
    incr state(depth)    

    if { $state(rootVisited) } {
	::tclu::abort "more than one rootelement; there can be only one !"
    }
    set state(rootVisited) 1
    set state(rootElem)    $name    

    lappend elemList $name
    lappend itemList $name


    $stackArr(currentElem) push $name
    #eval $code
    namespace eval schema $code
    $stackArr(currentElem) pop    

    incr state(depth) -1
    parseMsgOK_ $name
}

proc ::helpdoc::element {name code} {
    variable elemList
    variable itemList
    variable state
    variable stackArr
    variable elemArr

    parseMsg_ $name; puts ""
    incr state(depth)

    # check that $name does not exists
    if { [::tclu::lpresent $elemList $name] } {
	::tclu::abort "element \"$name\" already defined"
    }
    lappend elemList $name
    lappend itemList $name

    $stackArr(optional)   push 0
    $stackArr(interleave) push 0

    set parentElem [$stackArr(currentElem) peek]
    
    lappend elemArr(ELEMLIST,$parentElem)         $name
    lappend elemArr(OPTIONAL,$parentElem,$name)   [$stackArr(optional)   peek]
    lappend elemArr(INTERLEAVE,$parentElem,$name) [$stackArr(interleave) peek]
    lappend elemArr(REPETITION,$parentElem,$name) [$stackArr(repetition) peek]

    $stackArr(currentElem) push $name
    #eval $code
    namespace eval schema $code
    $stackArr(currentElem) pop

    $stackArr(optional)   pop
    $stackArr(interleave) pop

    incr state(depth) -1
    parseMsgOK_ $name
}

proc ::helpdoc::attribute {name code} { 
    # so far we assume attributes have arbitrary values (which means
    # we ignore code)

    variable itemList
    variable stackArr 
    variable elemArr

    parseMsg_ $name

    set currentElem [$stackArr(currentElem) peek]
    
    lappend itemList $name
    lappend elemArr(ATTRLIST,$currentElem) $name
    lappend attrArr(OPTIONAL,$currentElem) [$stackArr(optional) peek]
    
    puts ok
}


proc ::helpdoc::define {name code} {
    variable defineArr
    variable itemList

    parseMsg_ $name;
    
    lappend itemList $name
    set defineArr($name) $code

    puts ok
}


proc ::helpdoc::text {} { 
    # BEWARE: so far can be called only from element (because
    # attribute does not yet support ...)
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(TEXT,$currentElem) 1
}


proc ::helpdoc::String {} { 
    # BEWARE: so far can be called only from element (because
    # attribute does not yet support ...)
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(STRING,$currentElem) 1
}


proc ::helpdoc::clist {} { 
    # BEWARE: clist can be called only from element
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(CLIST,$currentElem) 1
}


proc ::helpdoc::ref {name} {    
    variable stackArr 
    variable elemArr
    variable defineArr

    parseMsg_ $name; 

    if { [info exists defineArr($name)] } {
	puts ""
	# the ref points to define, evaluate it
	#eval $defineArr($name)
	namespace eval schema $defineArr($name)	
	parseMsgOK_;	
	return
    }

    set currentElem [$stackArr(currentElem) peek]

    if { $currentElem != "" } {
	lappend elemArr(REFLIST,$currentElem) $name

	lappend elemArr(OPTIONAL,$currentElem,$name)   [$stackArr(optional)   peek]
	lappend elemArr(INTERLEAVE,$currentElem,$name) [$stackArr(interleave) peek]
	lappend elemArr(REPETITION,$currentElem,$name) [$stackArr(repetition) peek]
    } else {
	::tclu::abort "can't use \"ref\" outside element definition"
    }

    puts ok
}

proc ::helpdoc::ident {} {
    variable stackArr 
    variable elemArr
    
    set currentElem [$stackArr(currentElem) peek]

    if { $currentElem != "" } {
	set elemArr(IDENT,$currentElem) 1
    } else {
	::tclu::abort "can't use \"ident\" outside element definition"
    }  
}    

proc ::helpdoc::optional {code} {
    variable stackArr
    variable state

    parseMsg_; puts ""
    incr state(depth)

    $stackArr(optional) push 1
    # eval $code
    namespace eval schema $code
    $stackArr(optional) pop

    incr state(depth) -1
    parseMsgOK_ 
}

proc ::helpdoc::interleave {code} {
    variable stackArr
    variable state

    parseMsg_; puts "" 
    incr state(depth)

    $stackArr(interleave) push 1
    # eval $code
    namespace eval schema $code
    $stackArr(interleave) pop

    incr state(depth) -1
    parseMsgOK_ 
}

proc ::helpdoc::choice {code} {
    variable stackArr
    variable state

    # TODO: implement the CHOICE; so far this proc is dummy
    parseMsg_; puts ""
    incr state(depth)

    #eval $code
    namespace eval schema $code

    incr state(depth) -1
    parseMsgOK_ 
}

proc ::helpdoc::ancestorElements {} {
    parseMsg_ 
    # DO nothing (this means no validation for correctness will be done)
    puts ok
}

proc ::helpdoc::? {code} {
    repetition_ $code
}
proc ::helpdoc::* {code} {
    repetition_ $code
}
proc ::helpdoc::+ {code} {
    repetition_ $code
}
proc ::helpdoc::repetition_ {code} {
    variable stackArr
    variable state

    set type [tag -2]
    uplevel 1 "parseMsg_; puts {}"    
    incr state(depth)

    $stackArr(repetition) push $type
    #eval $code
    namespace eval schema $code
    $stackArr(repetition) pop    

    incr state(depth) -1
    uplevel 1 "parseMsgOK_"
}


proc ::helpdoc::assignRefs_ {} {
    variable elemList
    variable elemArr
 
    foreach elem $elemList {
	if { [info exists elemArr(REFLIST,$elem)] } {
	    # we have a ref
	    puts -nonewline "   $elem --> "
	    foreach ref $elemArr(REFLIST,$elem) {
		# check if ref points to "define"
		lappend elemArr(ELEMLIST,$elem) $ref

		puts -nonewline "$ref "
		
		# check that $ref exists
		if { ! [::tclu::lpresent $elemList $ref] } {
		    puts ""
		    ::tclu::abort "the \"$ref\" element has not been defined, yet it is referenced"
		}
	    }
	    puts ""
	}
    }
}


proc ::helpdoc::createTagCmds_ {} {
    variable state 
    variable elemList

    if { $state(rootElem) == {} } {
	::tclu::abort "rootelement was not defined"
    }

    # create the rootelement cmd
    puts "   creating $state(rootElem) cmd ... ok"
    proc ::helpdoc::tag::$state(rootElem) {args} {
	eval ::helpdoc::rootnameTag_ $args
    }

    # create all elements cmds

    foreach elem $elemList {
	if { $elem != $state(rootElem) } {
	    puts -nonewline "   creating $elem cmd ... "
	    proc ::helpdoc::tag::$elem {args} {
		eval ::helpdoc::elementTag_ $args
	    }
	    puts ok
	}
    }
}


# for the time being ...
proc helpdoc::parseMsg_ {{name {}}} {
    variable state
    
    set indent [::textutil::blank [expr (1+$state(depth)) * 3]]

    set tag [string toupper [tag -2]]
    puts -nonewline "${indent}parsing $tag $name ... "    
}

proc helpdoc::parseMsgOK_ {{name {}}} {
    variable state

    set indent [::textutil::blank [expr (1+$state(depth)) * 3]]

    set tag [string toupper [tag -2]]

    if { $name == "" } {
	puts "${indent}OK - parsing $tag completed"
    } else {
	puts "${indent}OK - parsing $tag $name completed"
    }
}
