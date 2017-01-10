#
# This file is sourced by "gen-emacs-mode" script and contains the
# machinery for generating the Emacs major mode files for Quantum
# ESPRESSO. It relies heavily on: (1) helpdoc utility, (2) INPUT_*.def
# files, and (3) templates provided in qe-emacs.templates/ directory.
#

###
# Purpose: extract the module name from INPUT_*.def filename
proc ::helpdoc::modulename_from_defname {deffile} {

    # strip the path from $deffile
    
    set def_file [file tail $deffile]
    
    # the module-name of INPUT_XX.def is [string tolower XX]

    return [string tolower [regsub {^INPUT_} [file rootname $def_file] {}]]
}

# Purpose: quote list, i.e., transform {a b c} to {"${prefix}a" "${prefix}b" "${prefix}c"},
#          where $prefix stands for the value of prefix
proc ::helpdoc::quote_list {lst {prefix {}}} {    
    set result ""
    foreach item $lst {
	append result "\"${prefix}$item\" "
    }
    return $result
}


# Purpose: store the namelists, variables, cards, etc. that are stated in the def file
#          (this proc is used by tree's walkproc method; see qe_mode_process_def)
proc ::helpdoc::getFontlockKeys {tree node action} {    
    variable arr
    variable inside_namelist    
    variable fontlock
    variable defun
    variable module

    set tag  [$tree get $node tag]
    set attr [getFromTree $tree $node attributes]

    catch {unset arr}
    attr2array_ arr $attr
    set name [arr name]

    if { $action == "enter" } {
	
	switch -- $tag {
	    supercard {
		set start [supercardStarttag]
		set end   [arr endtag]
		lappend defun($module,supercards) [list $start $end]
		
		if { $end ne {} } {
		    lappend fontlock(begin_supercards)   $start
		    lappend fontlock(end_supercards)     $end
		} else {
		    lappend fontlock(open_supercards)  $start
		}
	    }
	    namelist {
		set inside_namelist 1
		lappend fontlock(namelists) $name
		lappend defun($module,namelists) $name
	    }
	    var - dimension {
		if { $inside_namelist } {
		    set bare_name [lindex [split $name \(] 0]; # strips "(index)" from variable's name if var is a dimension
		    lappend fontlock(vars) $bare_name
		    lappend defun($module,vars) $bare_name
		    
		    # is the variable of string-type ?
		    set type [arr type]
		    if { "$type" eq "CHARACTER" } {
			#set fontlock(stringvar,$bare_name) 1
			set defun($module,stringvar,$bare_name) 1
		    }
		}
	    }
	    keyword {
		lappend fontlock(keywords) $name
		lappend defun($module,cards) $name; # make it easier (glue cards and keywords immediately)
	    }
	    card {	    
		set nameless [arr nameless]
		switch -- [string tolower $nameless] {
		    1 - true - yes - .true. {
			set name ""
		    }	    
		}
		if { $name != "" } {
		    lappend fontlock(cards) $name
		    lappend defun($module,cards) $name

		    set flags_txt [getDescendantText $tree $node flag enum]
		    set flags_use [getDescendantAttribute $tree $node flag use]
		    set flags     [regsub -all -- {\|} $flags_txt {}]
		    
		    if { $flags != "" } {
			append fontlock(flags) "$flags "
			#append defun($module,flags) "$flags "
			
			if { $flags_use == "optional" } {
			    #set fontlock(card_flags,$name) "\{ [string trim $flags_txt { }] \}"
			    set defun($module,card_flags,$name) "\{ [string trim $flags_txt { }] \}"
			} else {
			    #set fontlock(card_flags,$name) [string trim $flags_txt { }]
			    set defun($module,card_flags,$name) [string trim $flags_txt { }]
			}
		    }
		}
	    }
	}
    } else {
	if { $tag == "namelist" } {
	    set inside_namelist 0
	}
    }
}


# Purpose: initialize the machinery for generating qe-modes
proc ::helpdoc::qe_mode_init {} {
    variable init_qe_mode    
    global qe_modes_template_dir
    variable fontlock
    variable defun

    #catch {unset fontlock}
    #catch {unset defun}

    set init_qe_mode 1
    
    if { ! [info exists qe_modes_template_dir] } {
	# try with this
	puts stderr "
### variable \"qe_modes_template_dir\" is not defined ... aborting"
	exit 1
    }

    if { ! [file isdirectory $qe_modes_template_dir] } {
	puts stderr "
### variable \"qe_modes_template_dir\" points to nonexistent directory: $qe_modes_template_dir
"
	exit 1
    }

    # read the schema (and load tag's commands)

    readSchema
}


# Purpose: process a given *.def file for the pupose of generating qe-modes
proc ::helpdoc::qe_mode_process_def {deffile} {
    variable tree
    variable inside_namelist
    variable init_qe_mode

    if { ! [info exists init_qe_mode] || ! $init_qe_mode } {
	puts stderr "
### call ::helpdoc::qe_mode_init before the ::helpdoc::qe_mode_process_def
"
	exit 1
    }
    
    set inside_namelist 0
    
    # read the deffile

    puts "\n\n***\n*** Parsing definition file: $deffile\n***\n"
    namespace eval tag [list source $deffile]
    puts ""
    
    # walk through the $deffile and extract the info we need (i.e. names of namelists, variables, cards, and flags)
    
    $tree walkproc root -order both helpdoc::getFontlockKeys
}


# Purpose: generate requested qe-mode
proc ::helpdoc::qe_mode_generate {module_list} {
    variable fontlock
    variable defun
    variable opt
    global qe_modes_template_dir
    # load the variables needed by qe-mode.el.tcl
    
    set mode [string tolower $opt(mode)]
    
    if { $opt(modename) eq {} } {
	set modeName QE-$mode.x
    } else {
	set modeName $opt(modename)
    }

    # closed supercards, i.e., with starttag & endtag
    set sc_l [concat [value_of fontlock(begin_supercards)] [value_of fontlock(end_supercards)]]
    set closed_supercards [quote_list [lsort -unique $sc_l]]
    
    # open supercards, i.e., with only starttag
    set open_supercards  [quote_list [lsort -unique [value_of fontlock(open_supercards)]]]

    # cards & flags, namelists & vars
    set cards_l   [concat [value_of fontlock(keywords)] [value_of fontlock(cards)]]
    set cards     [quote_list [lsort -unique $cards_l]]
    set flags     [quote_list [lsort -unique [value_of fontlock(flags)]]]
    set namelists [quote_list [lsort -nocase -unique [value_of fontlock(namelists)]] $opt(nmlprefix)]
    set vars      [quote_list [lsort -nocase -unique [value_of fontlock(vars)]]]

    # load the templates
    
    set qe_template  [tclu::readFile [file join $qe_modes_template_dir qe-mode.el.tcl]]

    set header_template    [tclu::readFile [file join $qe_modes_template_dir header.el.tcl]] 
    set namelist_template  [tclu::readFile [file join $qe_modes_template_dir namelist.el.tcl]]
    set var_template       [tclu::readFile [file join $qe_modes_template_dir var.el.tcl]]
    set stringvar_template [tclu::readFile [file join $qe_modes_template_dir stringvar.el.tcl]]
    set card_template      [tclu::readFile [file join $qe_modes_template_dir card.el.tcl]]
    set card_noflags_template [tclu::readFile [file join $qe_modes_template_dir card-noflags.el.tcl]]

    set closed_supercard_template [tclu::readFile [file join $qe_modes_template_dir supercard.el.tcl]]
    set open_supercard_template   [tclu::readFile [file join $qe_modes_template_dir supercard-open.el.tcl]]
    
    # substitute the Tcl variables inside the qe-mode.el.tcl and write the resulting content

    set file $mode-mode.el
    puts "
writing Emacs major-mode file : $file
"
    set header [subst -nocommands -nobackslashes $header_template]    
    tclu::writeFile $file [subst -nocommands -nobackslashes $qe_template]
    
    # are we done ?
    
    if { ! $opt(funcs) } {
	# qe-mode function generation was not requested, we are done
	return
    }

    ########################################################################
    # make utility functions for each module in $Module_list

    set utility_functions ""

    puts "list-of-modules = $module_list\n"

    # create $module-insert-template functions for all modules for which the $module.in file exists

    foreach module $module_list {

	set in [file join $qe_modes_template_dir $module.in]
	
	if { [file exists $in] } {
	    set template [tclu::readFile $in]
	    set insert_template [tclu::readFile [file join $qe_modes_template_dir insert-template.el.tcl]]
	    
	    if { $utility_functions eq "" } {
		set utility_functions "
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; utility functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

"
	    }
	    append utility_functions [subst -nocommands -nobackslashes $insert_template]\n\n
	}
    }

    # make keyword functions for this $mode
    
    set keyword_functions ""

    foreach module $module_list {
	
	##################################################
	# create elisp functions for each namelist
	if { [value_of defun($module,namelists)] ne "" } {
	    
	    append keyword_functions "
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; $module- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

"
	    set module_nmls [lsort -nocase -unique [value_of defun($module,namelists)]]
	    
	    foreach nml $module_nmls {
		set namelist    "$opt(nmlprefix)[string toupper [string trimleft $nml &]]"
		set namelist_uc [string toupper [string trimleft $nml &]]
		append keyword_functions [subst -nocommands -nobackslashes $namelist_template]\n\n		
	    }
	}	

	##################################################
	# create elisp functions for each namelist's variable

	if { [value_of defun($module,vars)] ne "" } {
	    
	    append keyword_functions "
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; $module- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

"
	    set module_vars [lsort -nocase -unique [value_of defun($module,vars)]]

	    foreach var $module_vars {
	    
		set var_lc [string tolower $var]
	    
		# check if var is of string type
	    
		if { [info exists defun($module,stringvar,$var)] } {

		    set read_cmd read-string
		    
		    if { [string match {*dir} $var_lc] } {
			
			# variable probably want a directory name
			set read_cmd read-directory-name		
			
		    } elseif { [string match {*file*} $var_lc] } {
			
			# variable probably want a filename name
			set read_cmd read-file-name
		    } 
		    append keyword_functions [subst -nocommands -nobackslashes $stringvar_template]\n\n
		    
		} else {
		    
		    append keyword_functions [subst -nocommands -nobackslashes $var_template]\n\n
		}
	    }
	}

	##################################################
	# create elisp functions for each supercard
	
	if { [value_of defun($module,supercards)] ne "" } {
	    
	    append keyword_functions "	   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; $module- supercards functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

"
	    foreach start_end [lsort -unique $defun($module,supercards)] {
		
		set starttag [lindex $start_end 0]
		set endtag   [lindex $start_end 1]

		set supercard_uc [string toupper $starttag]

		if { $endtag ne {} } {
		    # it's a closed supercard		
		    append keyword_functions [subst -nocommands -nobackslashes $closed_supercard_template]\n\n
		} else {
		    # it's an open supercard		
		    append keyword_functions [subst -nocommands -nobackslashes $open_supercard_template]\n\n
		}
	    }
	}

	##################################################
	# create elisp functions for each card

	if { [value_of defun($module,cards)] ne "" } {

	    append keyword_functions "
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; $module- cards functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

"
	    set module_cards [lsort -unique $defun($module,cards)]
	    
	    foreach card $module_cards {
		
		set card_uc [string toupper $card]
		set card_lc [string tolower $card]
		
		if { [info exists defun($module,card_flags,$card)] } {
		    # it's a card with flags
		
		    set card_flags $defun($module,card_flags,$card)
		    append keyword_functions [subst -nocommands -nobackslashes $card_template]\n\n
		} else {
		    # it's a flagless card
		
		    append keyword_functions [subst -nocommands -nobackslashes $card_noflags_template]\n\n
		}
	    }
	}
    }
    
    set qe_func_template  [tclu::readFile [file join $qe_modes_template_dir qe-funcs.el.tcl]]

    set file qe-funcs.el
    puts "
writing function file : $file
"
    set header [subst -nocommands -nobackslashes $header_template]
    tclu::writeFile $file [subst -nocommands -nobackslashes $qe_func_template]
}


proc ::helpdoc::qe_master_generate {deffile_list} {
    variable opt
    global qe_modes_template_dir
    
    if { "$opt(mode)" ne "qe" } {
	puts stderr "
ERROR: so far master mode file can be generated only for mode = qe
"
	exit 1
    }

    puts "pwd = [pwd] 
"
    set header_template   [tclu::readFile [file join $qe_modes_template_dir header.el.tcl]] 
    set master_template   [tclu::readFile [file join $qe_modes_template_dir $opt(mode)-modes.el.tcl]]
    set specific_template [tclu::readFile [file join $qe_modes_template_dir autoload-specific.el.tcl]]

    set autoload_specific_modes ""

    foreach deffile $deffile_list {

	set module [modulename_from_defname $deffile]
	
	switch -exact -- $module {
	    atomic {
		set prog ld1.x
	    }
	    default {
		set prog $module.x
	    }
	}

	puts "mode file = $module-mode.el"
	
	if { [file exists $module-mode.el] } {
	    
	    # autolad only for those modules that have the specific mode defined
	    
	    append autoload_specific_modes [subst -nocommands -nobackslashes $specific_template]\n
	}
    }
    set file $opt(mode)-modes.el
    puts "

writing master mode file : $file
"
    set header [subst -nocommands -nobackslashes $header_template]
    tclu::writeFile $opt(mode)-modes.el [subst -nocommands -nobackslashes $master_template]
}
