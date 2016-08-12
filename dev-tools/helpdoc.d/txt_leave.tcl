variable indentNum

switch -exact -- $tag {
    
    input_description {	
    }
    
    intro {
	printf \n
    }
    
    toc {
    }
}

# simple elements

switch -exact -- $tag {

    info {
    }
    
    "default" {
    }
    
    status  {
    }
    
    label   {
    }
    
    message {
    }
    
    see     {
    }
    
    keyword {
    }
}


# composite elements

switch -exact -- $tag {
    list {
	if { [::tclu::lpresent $mode syntax] } {
	    syntaxFlush
	} else {
	    if { [printableVarDescription $tree $node] } {
		printf +--------------------------------------------------------------------
		printf \n
	    }
	}
    }
    
    format { # todo 
    }
    
    var - dimension - col - row {
	if { ! $vargroup && ! $dimensiongroup && ! $colgroup && ! $rowgroup && ! [::tclu::lpresent $mode syntax] } {
	    if { [printableVarDescription $tree $node] } {
		printf +--------------------------------------------------------------------\n
		set var_print 0
	    }
	}
    }
    
    vargroup - dimensiongroup - rowgroup - colgroup { # todo
	if { ! [::tclu::lpresent $mode syntax] } {
	    set $tag 0
	    if { [printableVarDescription $tree $node] } {
		printf +--------------------------------------------------------------------\n		
	    }
	}
    }
    
    table {
    }
    rows { 
	if { [::tclu::lpresent $mode syntax] && [::tclu::lpresent $mode rows] } {
	    manageRow 
	    manageRow 1	
	    manageRow 2
	    printRows	

	    unset rows
	    ::tclu::lpop mode
	}
    }
    cols { 
	if { [::tclu::lpresent $mode syntax] && [::tclu::lpresent $mode cols] } {
	    manageCol 
	    manageCol 1	
	    manageCol 2
	    printCols	

	    unset cols
	    ::tclu::lpop mode
	}
    }

    optional { 
	if { [::tclu::lpresent $mode rows] } {
	    append rows(line) "__optional::end__ "
	} elseif { [::tclu::lpresent $mode cols] } {
	    append cols(vline) "__optional::end__ "
	} elseif { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend "\}"
	}
    }
    conditional { 
	if { [::tclu::lpresent $mode rows] } {
	    append rows(line) "__conditional::end__ "
	} elseif { [::tclu::lpresent $mode cols] } {
	    append cols(vline) "__conditional::end__ "
	} elseif { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend "\]"
	}
    }
    group {
	incr txtDepth -1
	printf \\\\\\---\n
    }
    namelist {
	incr txtDepth -1
	printf "===END OF NAMELIST======================================================\n\n"
    }
    supercard {
	set name [supercardStarttag]
	set end  [arr endtag]

	if { $end   ne "" } { set name $name/$end }
	
	incr txtDepth -1
	set str "### END OF SUPERCARD :  $name "
	set l [string length $str]
	for {set i $l} {$i < 72} {incr i} {
	    append str #
	}
	
	printf $str\n\n
    }
    card {
    }
    linecard { 
	if { [::tclu::lpresent $mode syntax] } {
	    syntaxFlush
	}
    }
    flag { 
	if { [::tclu::lpresent $mode "description"] } {
	    printf +--------------------------------------------------------------------
	    puts $fid(txt) "\n" 
	}
    }
    enum { # todo
    }
    syntax {
	if { [::tclu::lpresent $mode syntax] } {	    	    
	    incr txtDepth -2
	    printf "\n/////////////////////////////////////////\n"
	    #printf "|______\n"
	    #printf "|                                       |"
	    #printf "+---------------------------------------+"
	    #printf "+---------------------------------------+\n"
	}
    }
    line { 
	if { [::tclu::lpresent $mode syntax] } {
	    syntaxFlush
	}
    }
    if {
	if { ! [::tclu::lpresent $mode description] } {	    
	    incr txtDepth -1
	    printf ENDIF
	}
    }
    choose { 
	if { ! [::tclu::lpresent $mode description] } {
	    printf ENDIF
	    printf ________________________________________________________________________\n
	}
    }
    when - elsewhen - otherwise {
	if { ! [::tclu::lpresent $mode description] } {
	    printf " "
	    incr txtDepth -1
	}
    }
    options {
	if { ! [::tclu::lpresent $mode syntax] } {
	    set options 0
	    set options_first 0
	}
    }
}


# some text structure stuff

switch -exact -- $tag {
    
    section {
	puts $fid(txt) ""
	incr txtDepth -1
    }
    subsection {
	puts $fid(txt) ""
	incr txtDepth -1
    }
    subsubsection {
	puts $fid(txt) ""
	incr txtDepth -1
    }
    paragraph {
    }
    text {
    }
}
