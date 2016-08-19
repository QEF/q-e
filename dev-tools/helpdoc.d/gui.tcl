#
# This file holds procs for checking the PWgui's modules against the
# INPUT_*.def files (and vice versa)
#


proc ::helpdoc::checkMsg {type msg} {
    puts [labelMsg ${type}: $msg]\n
}


proc ::helpdoc::checkGui_def_vs_module {} {
    variable def_item
    variable def_itemL
    variable moduleNewCode

    puts {
    -------------------------------------
     *** Checking DEF vs MODULE file ***
    -------------------------------------
    }
    
    foreach {name lowercase_name} $def_itemL {
	
	set def_type $def_item($name)	

	switch -- $def_type {
	    card {
		set def_mapping_type keyword
	    }
	    listvar - list {
		set def_mapping_type var
		set name [string trim $name ,]
	    }
	    default {
		set def_mapping_type ""
	    }
	}

	if { [info exists ::guib::moduleObj::module_item($name)] } {
	    
	    set module_type $::guib::moduleObj::module_item($name)
	    
	    if { $def_type != $module_type } {
		
		# take care of guib vs. helpdoc mappings
		
		if { $def_mapping_type != $module_type } {

		    set warning 1

		    # handle exceptions

		    switch -glob -- $name {
			first_image -
			intermediate_image -
			last_image {
			    if { $::module == "pw" } {
				# Don't report errors connected to: atomic_coordinates ...
				set warning 0
			    }
			}
		    }
		    if { $warning } {
			checkMsg WARNING "Type mismatch for item=$name.\n\tDef's type    = $def_type\n\tModule's type = $module_type"
		    }		    
		}
	    }
	    
	} else {
	    set module_name [module_getItemName_ $name]
	    if { $module_name != "" } {
		checkMsg WARNING "case-sensitivity mismatch for item $def_type $name.\n\tDef's name    = $name (type=$def_type)\n\tModule's name = $module_name (type=$module_type)"
	    } else {

		set error 1
		
		# handle exceptions
		
		switch -glob -- $name {
		    nwfts - test_wfs {
			if { $::module == "atomic" } {
			    # Don't report errors connected to: atomic_coordinates ...
			    set error 0
			}
		    }
		}
		if { $error } {
		    checkMsg ERROR "$def_type $name does not exist in MODULE file"

		    # let us attempt to privide default module's definition of variable

		    variable arr
		    attr2array_ arr $def_item($name,attrs)
		    set type [arr type]
		    
		    set options "   -label \"($name):\""
		    switch -glob -nocase -- $type {
			CHARACTER - STRING {
			    append options "
   -validate string"
			}
			LOGICAL {
			    append options "
   -widget    radiobox
   -textvalue { Yes No }	      
   -value     { .true. .false. }"
			}
			INT* {
			    append options "
   -validate int"
			}
			REAL {
			    append options "
   -validate fortranreal"
			}
		    }
		    
		    switch -exact -- $def_type {
			var - dimension {
			    append moduleNewCode "$def_type $name \{
$options
\}\n\n"
			}
		    }
		}
	    }
	}
    }
}

proc ::helpdoc::checkGui_module_vs_def {} {
    variable def_item
    variable def_itemL
    
    puts {
	
    -------------------------------------
     *** Checking MODULE vs DEF file ***
    -------------------------------------
    }
    
    foreach {name lowercase_name} $::guib::moduleObj::module_itemL {
	
	set module_type $::guib::moduleObj::module_item($name)	
	
	if { [info exists def_item($name)] } {
	    
	    set def_type $def_item($name)
	    	    
	    if { $def_type != $module_type } {
		
		# take care of guib vs. helpdoc mappings

		switch -- $def_type {
		    card {
			set def_mapping_type keyword
		    }
		    listvar {
			set def_mapping_type var
			set name [string trim $name ,]
		    }
		    default {
			set def_mapping_type ""
		    }
		}
		
		if { $def_mapping_type != $module_type } {

		    # handle exceptions

		    set warning 1
		    
		    switch -glob -- $name {
			first_image -
			intermediate_image -
			last_image {
			    if { $::module == "pw" } {
				# Don't report errors connected to: atomic_coordinates for pw.x ...
				set warning 0
			    }
			}
		    }

		    if { $warning } {
			checkMsg WARNING "Type mismatch for item=$name.\n\tModule's type = $module_type\n\tDef's type    = $def_type"
		    }
		}	    
	    }	    
	} else {
	    set def_name [def_getItemName $name]
	    if { $def_name != "" } {
		checkMsg WARNING "case-sensitivity mismatch for item $def_type $name.\n\tModule's name = $name (type=$module_type)\n\tDef's name    = $def_name (type=$def_type)"
	    } else {
		
		# handle exceptions
		
		set error 1

		switch -glob -- $name {
		    atomic_coordinates_* -
		    first_image -
		    intermediate_image -
		    last_image {
			if { $::module == "pw" } {
			    # Don't report errors connected to: atomic_coordinates ...
			    set error 0
			}
		    }
		    nwfts_* - test_wfs_* {
			if { $::module == "atomic" } {
			    set error 0
			}
		    }
		}
		if { $error }  {
		    checkMsg ERROR "$module_type $name does not exists in DEF file"
		}
	    }
	}
    }
}


#
# DEF's related proc's
#

proc ::helpdoc::def_loadDef {file} {
    variable tree
    variable def_item
    variable def_itemL

    if { [info exists def_item] } { unset def_item }
    if { [info exists def_itemL] } { unset def_itemL }
    
    # first read the schema (and load tag's commands)
    readSchema 
    
    # now read the file
    namespace eval tag [list source $file] 

    $tree walkproc root -order pre helpdoc::def_registerItems

    return $tree
}

proc ::helpdoc::def_checkExistance_ {tag name} {
    variable def_item

    set lowercase_name [string tolower $name]

    if { [info exists def_item(name,$lowercase_name)] } {
	puts [labelMsg WARNING: "item $name already exists (old-tag=$def_item(tag,$lowercase_name), new-tag=$tag).\nAutomatic checking is not reliable, please check item, $name, manually."]
    }
}

proc ::helpdoc::def_registerItem_ {tag name {attrs {}}} {
    variable def_item
    variable def_itemL

    def_checkExistance_ $tag $name
    
    set lowercase_name [string tolower $name]

    set def_item($name) $tag
    set def_item($name,attrs) $attrs
    
    append def_itemL "[def_addToItemList__ $name] "
}

proc ::helpdoc::def_addToItemList__ {name} {
    set lowercase_name [string tolower $name]
    return [list $name $lowercase_name]
}
proc ::helpdoc::def_getItemName {name} {
    variable def_itemL
    set lowercase_name [string tolower $name]
    foreach {Name LowercaseName} $def_itemL {
	if { $LowercaseName == $lowercase_name } {
	    return $Name
	}
    }
    return {}
}
proc ::helpdoc::def_getItemLowercaseName {name} {
    set lowercase_name [string tolower $name]
    foreach {Name LowercaseName} $def_itemL {
	if { $LowercaseName == $lowercase_name } {
	    return $lowercase_name
	}
    }
    return {}
}

proc ::helpdoc::def_registerItems {tree node action} {
    variable arr

    set tag  [$tree get $node tag]
    set attr [getFromTree $tree $node attributes]

    catch {unset arr}
    attr2array_ arr $attr
    set name [arr name]

    switch -- $tag {
	var - keyword - dimension - namelist - table {
	    def_registerItem_ $tag $name $attr
	}
	list {
	    def_registerItem_ $tag $name $attr

	    set names [getTextFromDescendant $tree $node format]
	    foreach name $names {
		def_registerItem_ listvar $name $attr
	    }
	}
	card {	    
	    set nameless [arr nameless]
	    switch -- [string tolower $nameless] {
		1 - true - yes - .true. { set name "" }	    
	    }
	    if { $name != "" } {
		def_registerItem_ $tag $name
	    }
	}
    }
}


#
# guib-MODULE's related procs
#

proc ::helpdoc::module_getItemName_ {name} {
    
    set lowercase_name [string tolower $name]
    
    foreach {Name LowercaseName} $::guib::moduleObj::module_itemL {
	if { $LowercaseName == $lowercase_name } {
	    return $Name
	}
    }
    return {}
}
