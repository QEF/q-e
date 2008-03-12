set dir [file dirname [info script]]
lappend auto_path $dir [file join $dir .. .. GUI Guib lib]

package require tclu 0.9
package require struct::tree  2.1
package require struct::stack 1.3
package require textutil

namespace eval ::helpdoc {
    variable dir [file dirname [info script]]

    # schema-related variables

    variable attrArr; # stores all about attributes
    variable elemArr; # stores all about elements
    variable defineArr; # stores all about define's

    variable elemList ""
    variable itemList ""

    variable state    
    array set state {
	depth       0
	rootVisited 0
	rootElem    ""
    }

    variable stackArr    
    array set stackArr [subst {
	repetition  [::struct::stack]
	optional    [::struct::stack]
	interleave  [::struct::stack]
	currentElem [::struct::stack]
    }]

    $stackArr(repetition)  push 1; # decimal-digit | + | * | ? (meaning integer-number of times, one-or-more, zero-or-more, zero-or-one)
    $stackArr(optional)    push 0
    $stackArr(interleave)  push 0
    $stackArr(currentElem) push ""

    # stack & tree for parsing input definitions

    variable tree  ""
    variable stack [::struct::stack]

    # output-related

    variable indentNum 3
    variable txtDepth 0
    variable fid
    variable head
    variable rbd_var 
    variable rbd_stack 
    variable rbd_info
    variable robodoc  [auto_execok robodoc]
    variable xsltproc [auto_execok xsltproc]
}

namespace eval ::helpdoc::tag {}
namespace eval ::helpdoc::schema {}
source [file join $::helpdoc::dir readSchema.tcl]

proc ::helpdoc::openOutputs {file} {
    variable fid 
    variable head

    set head [file rootname $file]

    set fid(xml) [open $head.xml w]
    # currently disabled formats
    #set fid(txt) [open $head.txt w]
    #set fid(rbd) [open $head.rbd w]
    
    puts $fid(xml) {<?xml version="1.0" encoding="ISO-8859-1"?>}
    puts $fid(xml) {<?xml-stylesheet type="text/xsl" href="input_xx.xsl"?>}
    puts $fid(xml) {<!-- FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST -->
    }

    #puts $fid(txt) "*** FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST ***\n"
    
    #puts $fid(rbd) "# *** FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST ***\n"
}

proc ::helpdoc::writeOutputs {} {
    variable tree 
    variable head 
    variable fid
    variable robodoc 
    variable xsltproc 
    variable rbd_info

    #$tree destroy
    
    puts ""
    foreach fmt [array names fid] {
	puts "File $head.$fmt has been written."
	close $fid($fmt)
    }
    
    # run XSLTPROC

    if { $xsltproc != "" } {
	catch [list exec $xsltproc $head.xml > $head.html]    
	puts "File $head.html has been written."
    }
    
    # run ROBODOC 

    if { 0 } {
	# currently disbabled

	if { $robodoc != "" } {
	    if { ! [file isdirectory $head.d] } {
		file mkdir $head.d
	    } else {
		foreach file [glob -nocomplain $head.d/*.html] {
		    file delete $file
		}
	    }
	    if { ! [file isdirectory $head.robodoc] } {
		file mkdir $head.robodoc
	    }
	    
	    file copy -force $head.rbd $head.robodoc/
	    
	    catch {exec $robodoc --doc $head.d/ --src $head.robodoc/ --documenttitle "Description of $rbd_info(program) input file"}
	    
	    if { [file exists $head.d/toc_index.html] } {
		file copy -force $head.d/toc_index.html $head.d/index.html
		puts "File $head.d/index.html has been written."
	    }	
	}
    }
}


proc ::helpdoc::readSchema {} {
    puts "\n***\n*** Parsing the helpdoc.schema\n***\n"
    namespace eval schema { ::source helpdoc.schema }

    puts "\n\n***\n*** Assigning ref's\n***\n"
    assignRefs_
    
    puts "\n\n***\n*** Creating tags commands\n***\n"
    createTagCmds_
}

proc ::helpdoc::process {fileList} {
    variable tree

    # first read the schema (and load tag's commands)
    readSchema 

    #puts "tag commands: [info procs ::helpdoc::tag::*]"

    foreach file $fileList {
	
	if { [file exists $file] } {	
	    
	    openOutputs $file	

	    puts "\n\n***\n*** Parsing definition file: $file\n***\n"
	    namespace eval tag [list source $file]	    
	    
	    $tree walkproc root -order both print
	    writeOutputs

	    $tree destroy
	} else {
	    puts stderr "file [file join [pwd] $file] does not exists : aborting ..."	    
	    exit 1
	}
    }
}
