namespace eval ::helpdoc::gui_help {
    variable helpContent
    variable helpNameList ""

    proc printHelp_ {channel} {
	variable helpContent
	variable helpNameList
	
	foreach name $helpNameList {
	    puts $channel "\n# ------------------------------------------------------------------------"	    
	    if { [llength $name] > 1 } {
		puts $channel "grouphelp [list $name] -helpfmt helpdoc -helptext [list $helpContent($name)]\n"
	    } else {
		puts $channel "help $name -helpfmt helpdoc -helptext [list $helpContent($name)]\n"
	    }
	}
    }

    proc addHelp_ {names helpTxt} {
	variable helpContent
	variable helpNameList
		
	::tclu::ladd helpNameList $names
	append helpContent($names) ${helpTxt}\n
    }

    proc if_exists_addHelp_ {name helpTxt} {
	if { [info exists ::guib::moduleObj::module_item($name)] } {
	    if { $::guib::moduleObj::module_item(ident,$name) != "" } {
		switch -- $::guib::moduleObj::module_item($name) {
		    var - dimension - table - text {
			# important: we must pass from name to ident		    
			
			addHelp_ $::guib::moduleObj::module_item(ident,$name) $helpTxt
			return 1
		    }
		}
	    }
	}
	return 0
    }

    proc grouphelp {names helpTxt} {
	foreach name $names {
	    
	    if { [info exists ::guib::moduleObj::module_item($name)] } {
		if { $::guib::moduleObj::module_item(ident,$name) != "" } {
		    switch -- $::guib::moduleObj::module_item($name) {
			var - dimension - table  {
			    lappend ok_names $::guib::moduleObj::module_item(ident,$name)
			}
			
		    }
		}
	    } else {
		# pehaps the variable name is a dimension-like name, i.e. dimen(i,j,k),
		# and there is the "(i,j,k)" mismatch between module and def name;
		# try to see if "dimen" only exists ...
		
		set name [lindex [split $name \(] 0]
		
		if { [info exists ::guib::moduleObj::module_item($name)] } {
		    if { $::guib::moduleObj::module_item(ident,$name) != "" } {
			switch -- $::guib::moduleObj::module_item($name) {
			    var - dimension - table  {
				lappend ok_names $::guib::moduleObj::module_item(ident,$name)
			    }
			    
			}
		    }
		}
	    }       	
	}
	if { [info exists ok_names] } {	    
	    addHelp_ $ok_names $helpTxt
	}
    }

    #proc grouphelp {names helpTxt} {
    #	foreach name $names {
    #	    if { ! [if_exists_addHelp_ $name $helpTxt] } {	     	
    #
    #		# pehaps the variable name is a dimension-like name, i.e. dimen(i,j,k),
    #		# and there is the "(i,j,k)" mismatch between module and def name;
    #		# try to see if "dimen" only exists ...
    #		puts "pre.name = $name"
    #		set name [lindex [split $name \(] 0]
    #		puts "post.name = $name"
    #		
    #		if_exists_addHelp_ $name $helpTxt
    #	    }
    #	}
    #}

    proc help {name helpTxt} {

	# hande exceptions
	
	switch -- $::module {
	    atomic {
		switch -- $name {		    
		    nwfts - test_wfs {
			# in module file we have nwfts_*
			#puts "[array names ::guib::moduleObj::module_item -glob ${name}_*]"
			set names [array names ::guib::moduleObj::module_item -glob ${name}_*]
			if { $names != "" } {
			    grouphelp $names $helpTxt			
			}
		    }
		}
	    }

	    ph {
		if { $name eq "alpha_mix(niter)" } {
		    # in module file we have alpha_mix(1)
		    set name alpha_mix(1)
		}

		# for the time being ...
		if { $name eq "dvscf_star" || $name eq "drho_star" } {
		    foreach elem {open dir ext basis pat} {
			lappend names $name%$elem
		    }
		    grouphelp $names $helpTxt
		    return
		}
	    }   		    
	}

	if { $name == "occupations_table" } {
	    puts "occupations_table"
	    puts "      def-exists [info exists ::helpdoc::def_item($name)]"
	    puts "   module-exists [info exists ::guib::moduleObj::module_item($name)]"
	    puts "    module-ident $::guib::moduleObj::module_item(ident,$name)"
	}

	if { ! [if_exists_addHelp_ $name $helpTxt] } {	     	

	    # pehaps the variable name is a dimension-like name, i.e. dimen(i,j,k),
	    # and there is the "(i,j,k)" mismatch between module and def name;
	    # try to see if "dimen" only exists ...

	    set name [lindex [split $name \(] 0]
	    
	    if_exists_addHelp_ $name $helpTxt
	}
    }
}


proc ::helpdoc::checkGui_makeHelpFile {deffile modulefile} {
    variable xsltproc
    variable helpfile
    variable xml_temp

    if { $xsltproc == "" } {
	::tclu::ERROR "can't find useable xsltproc, gui help file creation skipped"
    }
            
    # help file will be written to $helpfile

    set helpfile      [file tail [file rootname $modulefile]]-help.tcl
    set orig_helpfile [file rootname $modulefile]-help.tcl

    if { "$helpfile" == "$orig_helpfile" } {
	puts stderr [::tclu::labelMsg WARNING "file \"$orig_helpfile\" exists.\nMaking a $orig_helpfile.bak backup copy."]
	file copy -force $orig_helpfile $orig_helpfile.bak
    }


    # open/create a temporaty xml file ...
    
    set orig_xmlfile [file rootname $deffile].xml
    set xml_prefix   [file tail [file rootname $deffile]]

    if { "$xml_prefix.xml" == "$orig_xmlfile" } {
	# ups, we don't want to overwrite $xmlfile
	set xml_temp ${xml_prefix}_temp.xml
    } else {
	set xml_temp ${xml_prefix}.xml
    }
    set xml_fid [open $xml_temp w]


    # copy $orig_xmlfile to $xml_temp, but replace the stylesheet input_xx.xsl by guihelp.xsl

    ::tclu::lineread line $orig_xmlfile {
	if { [string match {<?xml-stylesheet*} $line] } {
	    puts $xml_fid {<?xml-stylesheet type="text/xsl" href="guihelp.xsl"?>}
	} else {
	    puts $xml_fid $line
	}
    }
    close $xml_fid
    puts "\n\tXml-file $xml_temp has been written.\n"

    catch [list exec $xsltproc $xml_temp > $xml_temp.tcl]    

    puts "\n\tAuxiliary help-file $xml_temp.tcl has been written.\n"


    # create a $helpfile

    namespace eval gui_help {

	set helpID [open $::helpdoc::helpfile w]

	puts $helpID {
#
# Help-file automatically created by helpdoc utility
#
#    !!! DO NOT EDIT: CHANGES WILL BE LOST !!!
#
	}
	
	source $::helpdoc::xml_temp.tcl
	printHelp_ $helpID
	close $helpID
    }
    puts "\n\tHelp-file $helpfile has been written.\n"
}

