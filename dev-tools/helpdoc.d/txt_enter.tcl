variable indentNum
set var_chars 15
set var_chars1 [expr $var_chars + 1]
set var_chars_indent [expr 15 + 5]

switch -exact -- $tag {
    
    input_description {
	if { [info exists ::opt(version)] && $::opt(version) ne {} } {
	    set arr(version) "(version: $::opt(version))"
	}
	printfNormalize [subst {
	    ------------------------------------------------------------------------
	    INPUT FILE DESCRIPTION
	    
	    Program: [arr program] / [arr package] / [arr distribution] [arr version]
	    ------------------------------------------------------------------------
	}]
	printf \n
    }
    
    intro {
	printf $content\n
    }
    
    toc {
	# o-la-la ...
    }
}

# simple elements

switch -exact -- $tag {
    
    label   {
	if { ! [::tclu::lpresent $mode description] } {
	    printf [string toupper $content]\n
	}
    }
    
    message {
	if { ! [::tclu::lpresent $mode description] } {
	    printf $content\n	    
	}
    }
        
    keyword {
	if { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend [arr name]
	}
    }
}


if { ! $vargroup && ! $dimensiongroup && ! $colgroup && ! $rowgroup && ! [::tclu::lpresent $mode syntax] } {

    switch -exact -- $tag {
	
	info {
	    if { $options == 0 } {
		printf [labelMsg [format "%-${var_chars}s" Description:] $content]
	    } else {
		if { ! $options_first } {		    
		    printf " "
		}
		printf [labelMsg [format "%-${var_chars}s" {}] $content]
		set options_first 0
	    }
	}    

	opt {
	    if { [string trim $content] ne {} } {
		if { ! $options_first } {
		    printf " "
		}
		printf [labelMsg [format "%-${var_chars}s" {}] "[arr val] :"]
		printf [labelMsg [format "%-${var_chars_indent}s" {}] $content]
	    } else {
		printf [labelMsg [format "%-${var_chars_indent}s" {}] "[arr val]"]
	    }
	    set options_first 0
	}
	
	"default" {
	    printf [labelMsg [format "%-${var_chars}s" Default:] $content]
	}
    
	status  {
	    printf [labelMsg [format "%-${var_chars}s" Status:] $content]
	}

	see     {
	    printf [labelMsg [format "%-${var_chars}s" See:] $content]
	}	
    }
}


# composite elements


switch -exact -- $tag {
    options {
	if { ! [::tclu::lpresent $mode syntax] } {	
	    set options 1
	    set options_first 1
	    printf [format "%-${var_chars}s" Description:]
	}
    }

    var - col - row {
	if { ! $vargroup && ! $colgroup && ! $rowgroup && ! [::tclu::lpresent $mode syntax] } {	    

	    if { [printableVarDescription $tree $node] } {
		printf +--------------------------------------------------------------------
		printf [labelMsg [format "%-${var_chars}s" Variable:] [arr name]]\n
		printf [labelMsg [format "%-${var_chars}s" Type:] [arr type]]	    
	    }
	}

	if { $tag == "var" && [::tclu::lpresent $mode syntax] } {
	    syntaxAppend [arr name]
	}
    }
    
    dimension {
	if { ! $dimensiongroup && ! [::tclu::lpresent $mode syntax] } {
	    if { [printableVarDescription $tree $node] } {
		printf +--------------------------------------------------------------------
		printf [labelMsg [format "%-${var_chars}s" Variable:] "[arr name](i), i=[arr start],[arr end]"]\n
		printf [labelMsg [format "%-${var_chars}s" Type:] [arr type]]
	    }
	}

	if { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend "[arr name], i=[arr start],[arr end]"
	}
    }

    vargroup - dimensiongroup - colgroup - rowgroup {
	if { ($tag == "colgroup" || $tag == "rowgroup") && ! [::tclu::lpresent $mode description] } {
	    return
	}
	if { ! [::tclu::lpresent $mode syntax] } {
	    set $tag 1
	    foreach child [$tree descendants $node] {
		set _tag  [getFromTree $tree $child tag]
		set _attr [getFromTree $tree $child attributes]
		set _text [getFromTree $tree $child text]
		
		attr2array_ _arr $_attr
		
		switch -exact -- $_tag {
		    var - col - row {
			append Data(vars) "$_arr(name), "
		    }
		    dimension {
			append Data(dims) "${_arr(name)}(i), "
		    }		
		    status - "default" - info - see {
			set Data($_tag) [formatString $_text]
		    }
		}
	    }
	    
	    if { [printableVarDescription $tree $node] } {
		printf +--------------------------------------------------------------------
		if { $tag != "dimensiongroup" } {
		    printf [labelMsg [format "%-${var_chars}s" Variables:] [string trim $Data(vars) {, }]]\n
		} else {
		    printf [labelMsg [format "%-${var_chars}s" Variables:] "${Data(dims)}i=[arr start],[arr end]"]\n	    
		}
		printf [labelMsg [format "%-${var_chars}s" Type:] [arr type]]
		foreach field {default status see info} {	    
		    if { [info exists Data($field)] } {
			if { $field != "info" } {
			    set label [string totitle $field]:
			} else {
			    set label Description:
			}
			printf [labelMsg [format "%-${var_chars}s" $label] $Data($field)]
		    }
		}
	    }
	}        
    }
}

switch -exact -- $tag {
    list { 
	if { ! [::tclu::lpresent $mode syntax] } {
	    if { [printableVarDescription $tree $node] } {
		set vars [getTextFromDescendant $tree $node format]
		printf +--------------------------------------------------------------------
		printf [labelMsg [format "%-${var_chars}s" Variables:] $vars]\n
		printf [labelMsg [format "%-${var_chars}s" Type:] [arr type]]
	    }
	}

	if { $tag == "var" && [::tclu::lpresent $mode syntax] } {
	    syntaxAppend [arr name]
	}

    }
    
    format { 
	if { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend $content
	}
    }    
    table {
    }
    rows { 
	if { [::tclu::lpresent $mode syntax] } {	    
	    set rows(start) [arr start]
	    set rows(end)   [arr end]
	    lappend mode rows
	}
    }
    cols {
	if { [::tclu::lpresent $mode syntax] } {	    
	    set cols(start) [arr start]
	    set cols(end)   [arr end]
	    lappend mode cols
	}
    }
    rowgroup { 
	set rowgroup 1
    }
    colgroup { 
	set colgroup 1
    }
    col { 
	if { [::tclu::lpresent $mode rows] } {
	    append rows(line) "[arr name] "
	}
    }
    row { 
	if { [::tclu::lpresent $mode cols] } {
	    append cols(vline) "[arr name] "
	}
    }
    optional    { 
	if { [::tclu::lpresent $mode rows] } {
	    append rows(line) "__optional::begin__ "
	} elseif { [::tclu::lpresent $mode cols] } {
	    append cols(vline) "__optional::begin__ "
	} elseif { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend "\{"
	}
    }
    conditional {
	if { [::tclu::lpresent $mode rows] } {
	    append rows(line) "__conditional::begin__ "
	} elseif { [::tclu::lpresent $mode cols] } {
	    append cols(vline) "__conditional::begin__ "
	} elseif { [::tclu::lpresent $mode syntax] } {
	    syntaxAppend "\["
	}
    }
    group { # todo
	printf ///---
	incr txtDepth
    }
    namelist {
	printf ========================================================================
	printf "NAMELIST: &[arr name]\n"
	incr txtDepth
    }
    supercard {
	# rule is as follows: if -startag and -endtag exists use then, otherwise use the name instead.
	set start [supercardStarttag]
	set name  $start
	set end   [arr endtag]
	set rem   [formatString [arr remark]]
	
	if { $end ne "" } { set name $name/$end }

	printf "########################################################################"
	printf "| SUPERCARD: $name"
	if { $end != "" } {
	    printf "| this supercard is enclosed within the keywords:"
	    printf "|"
	    printf "| $start"
	    printf "|    ... content of the supercard here ..."
	    printf "| $end"
	    
	} else {
	    printf "| this supercard starts with the keyword:"
	    printf "|"
	    printf "| $start"
	    printf "|    ... content of the supercard here ..."
	}
	printf "|"
	if { $rem ne "" } {
	    printf "| REMARK:"
	    foreach line [split $rem \n] {
		printf "| $rem"
	    }
	    printf "|"
	}
	printf "| The syntax of supercard's content follows below:\n"
	incr txtDepth
    }
    card {
	if { ! [::tclu::lpresent $mode card] } {

	    lappend mode card

	    set flags [string trim [getTextFromDescendant $tree $node enum]]
	    set use   [getAttributeFromDescendantPath $tree $node flag use]
	    
	    if { $use == "optional" } {
		set flag "{ $flags }"
	    } else {
		set flag "$flags"
	    } 
	    
	    set card(name) [arr name]
	    set card(flag) $flag

	    set nameless [arr nameless]
	    switch -- [string tolower $nameless] {
		1 - true - yes - .true. { set card(name) "" }	    
	    }
	    
	    printf ========================================================================
	    printf "CARD: $card(name) $flag\n"
	    incr txtDepth
	    
	    # first parse subtree in syntax mode

	    txt_subtree $tree $node syntax
	    
	    # now parse subtree in description mode

	    # card's items will be described only when at least one of
	    # info, status or see records is present.
	    
	    set Info   [getTextFromDescendant $tree $node info]
	    set Status [getTextFromDescendant $tree $node status]
	    set See    [getTextFromDescendant $tree $node see]
	    set Opt    [getTextFromDescendant $tree $node opt]

	    if { $Info != "" || $Status != "" || $See != "" || $Opt != "" } {
		printf "DESCRIPTION OF ITEMS:\n"
		incr txtDepth	    
		txt_subtree $tree $node description
		incr txtDepth -1
	    }

	    incr txtDepth -1
	    printf "===END OF CARD==========================================================\n\n"
	  
	    ::tclu::lpop mode
	    ::struct::tree::prune	    
	}
    }
    linecard { 
	if { ! [::tclu::lpresent $mode card] } {

	    lappend mode card
	    
	    set card(name) ""
	    set card(flag) ""

	    printf ========================================================================
	    printf "Line of input:\n"
	    incr txtDepth
	    
	    # first parse subtree in syntax mode

	    incr txtDepth
	    txt_subtree $tree $node syntax
	    incr txtDepth -1
	    printf \n
	    
	    # now parse subtree in description mode
	    
	    printf "DESCRIPTION OF ITEMS:\n"
	    incr txtDepth
	    
	    txt_subtree $tree $node description
	    
	    incr txtDepth -2
	    printf "===End of line-of-input=================================================\n\n"
	    
	    ::tclu::lpop mode
	    ::struct::tree::prune	    
	}
    }
    flag { 
	if { ! [::tclu::lpresent $mode syntax] } {
	    printf +--------------------------------------------------------------------
	    printf [labelMsg [format "%-${var_chars}s" "Card's flags:"] $card(flag)]\n
	}
    }
    enum { 
	# nothing
    }
    syntax { 
	if { [::tclu::lpresent $mode syntax] } {	    
	    set _flags [arr flag]
	    if { $_flags == "" } {
		set flags $card(flag)
	    } else {
		set flags $_flags
	    }
	    
	    printf "/////////////////////////////////////////"
	    printf "// Syntax:                             //"
	    printf "/////////////////////////////////////////\n"

	    incr txtDepth
	    if { $card(name) != "" } {
		printf "$card(name) $flags"
	    }
	    incr txtDepth	    
	}
    }
    line { 
	# nothing ??
    }
    if { 
	if { ! [::tclu::lpresent $mode description] } {
	    printf "* IF [arr test] : \n"
	    incr txtDepth
	}
    }
    choose {
	if { ! [::tclu::lpresent $mode description] } {
	    printf ________________________________________________________________________
	}
    }
    when { 
	if { ! [::tclu::lpresent $mode description] } {
	    printf "* IF [arr test] : \n"
	    incr txtDepth
	}
    }
    elsewhen { 
	if { ! [::tclu::lpresent $mode description] } {
	    printf "* ELSE IF [arr test] : \n"
	    incr txtDepth
	}
    }
    otherwise {
	if { ! [::tclu::lpresent $mode description] } {
	    printf "* ELSE : \n"
	    incr txtDepth
	}
    }
}


# some text structure stuff

switch -exact -- $tag {
    
    section {
	printf "\n:::: [arr title]\n"
	incr txtDepth
    }
    subsection {
	printf "\n::: [arr title]\n"
	incr txtDepth
    }
    subsubsection {
	printf "\n:: [arr title]\n"
	incr txtDepth
    }
    paragraph {
	printf "* [arr title]\n"
    }
    text {
	printf $content\n
    }
}
