set dir [file dirname [info script]]

# do we use the svn version

set libdir [file join $dir .. .. GUI Guib lib]
#
if { [file isdirectory $libdir] } {    
    lappend auto_path $dir $libdir

} else {
    # do we use the tarball version

    set libdir [glob -nocomplain -directory [file join $dir .. ..]  \
                    [file join PWgui-* lib Guib-* lib]]
    #
    if { $libdir != {} && [file isdirectory $libdir] } {
        lappend auto_path $dir $libdir

    } else {
        puts stderr "no usable directory containing the tclu library found; you need to download the PWgui package"
        exit 1
    }
}

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

    # TXT variables
    variable vargroup 0
    variable dimensiongroup 0
    variable colgroup 0
    variable rowgroup 0
    variable options 0
    variable options_first 0
}

namespace eval ::helpdoc::tag {}
namespace eval ::helpdoc::schema {}
source [file join $::helpdoc::dir readSchema.tcl]

proc ::helpdoc::openOutputs {file} {
    variable fid 
    variable head

    set head [file rootname $file]

    set fid(xml) [open $head.xml w]
    set fid(txt) [open $head.txt w]

    # currently disabled
    #set fid(rbd) [open $head.rbd w]
    
    puts $fid(xml) {<?xml version="1.0" encoding="ISO-8859-1"?>}
    puts $fid(xml) {<?xml-stylesheet type="text/xsl" href="input_xx.xsl"?>}
    puts $fid(xml) {<!-- FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST -->
    }

    puts $fid(txt) "*** FILE AUTOMATICALLY CREATED: DO NOT EDIT, CHANGES WILL BE LOST ***\n"
    
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

    # TXT
    puts $fid(txt) "This file has been created by helpdoc utility on [clock format [clock seconds]]"

    puts ""
    foreach fmt [array names fid] {
	puts "File $head.$fmt has been written."
	close $fid($fmt)
    }
    flush stdout

    # run XSLTPROC

    if { $xsltproc != "" } {
	puts -nonewline "   Executing:  $xsltproc --stringparam version \"$::opt(version)\" --stringparam current-date \"[clock format [clock seconds]]\" $head.xml > $head.html ..."
	
	if { [catch [list exec $xsltproc --stringparam version "$::opt(version)" --stringparam current-date [clock format [clock seconds]] $head.xml > $head.html] errorMsg] } {
	    puts " \[Error\]"
	    puts "Execution of xsltproc failed with error message:\n\n$errorMsg"
	} else {
	    puts " \[OK\]"
	    puts "File $head.html has been written."
	}	
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
    namespace eval schema { ::source [file join $basedir helpdoc.schema] }

    puts "\n\n***\n*** Assigning ref's\n***\n"
    assignRefs_
    
    puts "\n\n***\n*** Creating tags commands\n***\n"
    createTagCmds_

    #puts "\n\n***\n*** Creating regexps for @-tags\n***\n"
    createAtcmdRegexp_
}



proc ::helpdoc::print_xml {tree node action} {
    variable fid
    
    set depth [$tree depth $node]

    set tag        [$tree get $node tag]
    set attributes [getFromTree $tree $node attributes]
    set content    [getFromTree $tree $node text]
       
    xml_tag_${action} $tag $attributes $content $depth
}

proc ::helpdoc::print_txt {tree node action} {
    variable fid
    
    set depth [$tree depth $node]

    set tag        [$tree get $node tag]
    set attributes [getFromTree $tree $node attributes]
    set content    [getFromTree $tree $node text]

    

    txt_tag_${action} $tree $node $tag $attributes $content [expr $depth - 1]

    # currently disabled:

    # robodoc
    #rbd_tag_${action} $tag $attributes $content $depth
}



proc ::helpdoc::process {fileList} {
    variable tree
    variable vargroup 
    variable dimensiongroup 
    variable mode

    # first read the schema (and load tag's commands)
    readSchema 

    #puts "tag commands: [info procs ::helpdoc::tag::*]"

    foreach file $fileList {

	set vargroup 0
	set dimensiongroup 0
	
	if { [file exists $file] } {	
	    
	    openOutputs $file	

	    puts "\n\n***\n*** Parsing definition file: $file\n***\n"
	    namespace eval tag [list source $file]	    
	    
	    set mode default

	    $tree walkproc root -order both print_xml
	    $tree walkproc root -order both print_txt
	    writeOutputs

	    $tree destroy
	    unset mode
	} else {
	    puts stderr "file [file join [pwd] $file] does not exists : aborting ..."	    
	    exit 1
	}
    }
}
