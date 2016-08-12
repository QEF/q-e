
namespace eval ::helpdoc::schema {

    # here is the definition of Tcl-commands that are used in schema

    proc rootelement {name code} { uplevel 1 [list ::helpdoc::rootelement $name $code] }
    proc element     {name code} { uplevel 1 [list ::helpdoc::element $name $code] }
    proc attribute   {name code} { uplevel 1 [list ::helpdoc::attribute $name $code] }
    proc @           {name code} { uplevel 1 [list ::helpdoc::@ $name $code] }
    
    proc word        {}          { uplevel 1 [list ::helpdoc::word] }
    proc string      {}          { uplevel 1 [list ::helpdoc::String] }
    proc text        {}          { uplevel 1 [list ::helpdoc::text] }
    proc clist       {}          { uplevel 1 [list ::helpdoc::clist] }
    proc plist       {}          { uplevel 1 [list ::helpdoc::plist] }

    proc ident       {}          { uplevel 1 [list ::helpdoc::ident] }
    
    proc define      {name code} { uplevel 1 [list ::helpdoc::define $name $code] }    
    proc ref         {name}      { uplevel 1 [list ::helpdoc::ref $name] }

    proc optional    {code}      { uplevel 1 [list ::helpdoc::optional $code] }    
    proc interleave  {code}      { uplevel 1 [list ::helpdoc::interleave $code] }
    proc choice      {code}      { uplevel 1 [list ::helpdoc::choice $code] }
    proc ?           {code}      { uplevel 1 [list ::helpdoc::? $code] }
    proc *           {code}      { uplevel 1 [list ::helpdoc::* $code] }
    proc +           {code}      { uplevel 1 [list ::helpdoc::+ $code] }

    proc ancestorElements {}     { uplevel 1 [list ::helpdoc::ancestorElements] }
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

proc ::helpdoc::@ {name code} {
    # PURPOSE: used to define the @-prefixed tags
    variable atTags
    
    parseMsg_ $name;

    # let's register the @ command at its type ...

    set type  [lindex $code 0]

    switch -exact $type {
	empty {
	    lappend atTags(empty) $name	    
	}
	special {
	    set regex [lindex $code 1]
	    if { $regex eq {} } {
		::tclu::abort "no regexp specified for special type @-prefixed command @name"
	    }
	    lappend atTags(special) $name
	    set atTags(regexp,$name) $regex
	}
	word {
	    lappend atTags(word) $name
	}
	varname {
	    lappend atTags(varname) $name
	}
	text - list - string - clist - plist {
	    lappend atTags(text) $name
	}
	default {
	    ::tclu::abort "unsupported content/datatype ($code) for @-prefixed command $name, must be obe of empty, word, varname, text, or special"
	}
    }
    
    puts ok;
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


proc ::helpdoc::empty {} {    
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(EMPTY,$currentElem) 1; # so far not used ...
}

proc ::helpdoc::word {} { 
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(WORD,$currentElem) 1
}

proc ::helpdoc::String {} { 
    # BEWARE: so far can be called only from element (because
    # attribute does not yet support ...)
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(STRING,$currentElem) 1
}


proc ::helpdoc::text {} { 
    # BEWARE: so far can be called only from element (because
    # attribute does not yet support ...)
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(TEXT,$currentElem) 1
}


proc ::helpdoc::clist {} { 
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(CLIST,$currentElem) 1
}


proc ::helpdoc::plist {} { 
    variable stackArr 
    variable elemArr
    set currentElem [$stackArr(currentElem) peek]
    set elemArr(PLIST,$currentElem) 1
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


proc ::helpdoc::createAtcmdRegexp_ {} {
    # PURPOSE: create regexps that will be used to properly substitute all @tags within correponding XML or TXT analogues ...
    variable atTags
    variable xml_re
    
    if { ! [array exists atTags] } {
	return
    }

    # create regexps
    
    # empty 
    if { [info exists atTags(empty)] } {

	set tag_re [join $atTags(empty) |]
	set xml_re(empty) "@($tag_re)(?=\\s+)"
	#puts "empty tags regexp: $xml_re(empty)"
    }

    # word
    if { [info exists atTags(word)] } {

	set tag_re [join $atTags(word) |]
	append xml_re(word) "@($tag_re)"
	#                     1
	append xml_re(word) {\s+(\S*[^\s\.\?\!\;\),]+)}; # ".!?,:;)" are excluded as last characters from word
	#                       2
	#puts "word tags regexp: $xml_re(word)"
    }
    
    # varname
    if { [info exists atTags(varname)] } {

	set tag_re [join $atTags(varname) |]
	append xml_re(varname) "@($tag_re)"
	#                        1
	append xml_re(varname) {\s+(\w+([%]\w+)*)}
	#                          2   3
	#puts "varname tags regexp: $xml_re(varname)"
    }

    # text
    if { [info exists atTags(text)] } {
	
	set tag_re [join $atTags(text) |]

	set xml_re(text,tagsonly) "@($tag_re)"
	append xml_re(text) "@($tag_re)"	
	#                     1
	append xml_re(text) {\s+((?:\{)((?:\{.*\}|[^\{])*)(?:\})|([^\s<]+))}
	#                       2      3                         4
	#puts "text tags regexp: $xml_re(text)"
    }
}


proc ::helpdoc::xml_atTags {content} {
    # PURPOSE: substitute all instances of @tag within the $content with the XML's <tag></tag> or <tag/>
    variable xml_re

    # empty
    if { [info exists xml_re(empty)] } {
	
	set content [regsub -all $xml_re(empty) $content {<\1/>}]
    }
    
    # varname & word
    foreach item {varname word} {
	
	if { [info exists xml_re($item)] } {	    
	    set content [regsub -all "$xml_re($item)" $content {<\1>\2</\1>}]
	}
    }

    # text
    if { [info exists xml_re(text)] } {

	while { [regexp $xml_re(text,tagsonly) $content] } {
	    set content [regsub -all $xml_re(text) $content {<\1>\3\4</\1>}]

	    # safety check: prevent infinite loop if something goes wrong ...
	    incr i
	    if { $i > 100 } {break}
	}
    }
    return $content
}


proc ::helpdoc::txt_atTags {content} {
    # PURPOSE: either ignore all specially treat all instances of
    # @tag's within the $content with as to get read of @tags in the
    # generated TXT representation
    variable xml_re

    # special processing for tag: hr

    set content [regsub -all {(@hr)} $content \
		     _____________________________________________________________________]
    

    # other tags ...
    
    # empty
    if { [info exists xml_re(empty)] } {
	
	set content [regsub -all $xml_re(empty) $content {}]
    }
    
    # varname & word
    foreach item {varname word} {
	
	if { [info exists xml_re($item)] } {	    
	    set content [regsub -all "$xml_re($item)" $content {\2}]
	}
    }
    
    # text
    if { [info exists xml_re(text)] } {
	
	while { [regexp $xml_re(text,tagsonly) $content] } {
	    set content [regsub -all $xml_re(text) $content {\3\4}]
	    
	    # safety check: prevent infinite loop if something goes wrong ...
	    incr i
	    if { $i > 100 } {break}
	}
    }
    
    return $content
}
