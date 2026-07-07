#
# TXT
#

proc ::helpdoc::attr2array_ {arrayName attributes} {
    upvar $arrayName arrayVar

    catch {array unset arrayVar}; # EXPERIMENTAL: this should be the
			          # desired behavior, because one wants to
			          # tranform attribute-list to associative
			          # array, hence the previous key-value
			          # pairs should be cleared
    
    foreach {name value} [::textutil::splitx $attributes "=\"|\"\[ \n\r\\t\]|\"$"] {
	if { $name != "" } {
	    set arrayVar($name) [string trim $value =]
	}
    }
}

proc ::helpdoc::arr {elem} {
    variable arr

    if { [info exists arr($elem)] } {
	return $arr($elem)
    } 
    return ""
}

proc ::helpdoc::printf {content {extraSpace 0}} {
    variable txtDepth
    variable indentNum
    variable fid

    set indent [indent $txtDepth]
    if { $extraSpace > 0 } {
	set indent $indent[::textutil::blank $extraSpace]
    }
    foreach line [split $content \n] {
	puts $fid(txt) ${indent}$line
    }
}

proc helpdoc::printfNormalize {content} {
    variable txtDepth
    variable fid
    
    puts $fid(txt) [formatString $content $txtDepth]
}


proc helpdoc::labelMsg {label msg} {
    set il 1
    set len [string length $label]
    set message {}
    foreach line [split [string trim $msg] \n] {
        if { $il == 1 } {
            append message [::format "%${len}s %s" $label $line]
            incr il
        } else {
            append message [::format "\n%${len}s %s" {} $line]
        }
    }
    return $message
}

proc ::helpdoc::unitsNotation_ {expr} {
    # Render a units product-of-powers expression in conventional notation,
    # e.g. "Ry bohr^-1" -> "Ry/bohr". Keep in sync with the "units-notation"
    # template in dev-tools/input_xx.xsl.
    set num {}
    set den {}
    foreach tok [split [string trim $expr]] {
        if { $tok eq {} } { continue }
        if { [regexp {^(.+)\^-(.+)$} $tok -> base exp] } {
            if { $exp eq "1" } {
                lappend den $base
            } else {
                lappend den "$base^$exp"
            }
        } else {
            lappend num $tok
        }
    }
    set out [join $num " "]
    foreach d $den {
        append out "/$d"
    }
    return $out
}

proc ::helpdoc::unitsGloss {expr {kind dim}} {
    # Render a <units> or <dimensionality> value as a human-readable phrase.
    # kind=dim names the physical quantity; kind=units renders the unit in
    # conventional notation (atomic-unit composites name the unit system).
    # Keep the lookup tables in sync with the "units-gloss" template in
    # dev-tools/input_xx.xsl.
    set u [string trim $expr]
    if { $u eq {} } { return $u }
    set key [join [split $u] " "]

    if { $kind eq "units" } {
        # atomic-unit-system composites: name the unit system only
        array set sys {
            {bohr electron_mass^1/2 Ry^-1/2}      {Rydberg atomic units}
            {bohr electron_mass^1/2 Hartree^-1/2} {Hartree atomic units}
            {Ry e^-1 bohr^-1}                     {Rydberg atomic units}
            {Hartree e^-1 bohr^-1}                {Hartree atomic units}
        }
        if { [info exists sys($key)] } {
            return $sys($key)
        }
        # units: render the unit in conventional notation, no quantity name
        return [unitsNotation_ $key]
    }

    # dimensionality: name the physical quantity
    array set gloss {
        {Ry bohr^-1}                          {force (Ry/bohr)}
        {Hartree bohr^-1}                     {force (Hartree/bohr)}
        {Ry bohr^-3}                          {pressure (Ry/bohr^3)}
        {states eV^-1}                        {states/eV}
        {energy length^-1}                    {force}
        {energy length^-3}                    {pressure}
        {time^-1}                             {frequency}
        {length time^-1}                      {velocity}
        {energy charge^-1}                    {electric potential}
        {energy charge^-1 length^-1}          {electric field}
        {charge length^-3}                    {charge density}
    }
    if { [info exists gloss($key)] } {
        return $gloss($key)
    }
    # single token or unrecognized composite: fall back to the raw expression
    return $u
}

proc ::helpdoc::unitsGlossClist {content {kind dim}} {
    # Render the content of a <units>/<dimensionality> element: a single
    # product-of-powers expression, or a keyed comma-list whose keyed entries
    # are glossed to "<value> for index <i>" / "... otherwise".
    set c [string trim $content]
    if { $c eq {} } { return $c }
    if { [string first , $c] < 0 && [string first : $c] < 0 } {
        # not a comma-list and not keyed: a single plain expression
        return [unitsGloss $c $kind]
    }
    set keyed {}
    set unkeyed {}
    foreach tok [split $c ,] {
        set tok [string trim $tok]
        if { $tok eq {} } { continue }
        if { [regexp {^([0-9]+(?:-[0-9]+)?)\s*:\s*(.*)$} $tok -> key val] } {
            lappend keyed [unitsGlossKeyed_ $key [unitsGloss [string trim $val] $kind]]
        } else {
            lappend unkeyed [unitsGloss $tok $kind]
        }
    }
    set out $keyed
    foreach u $unkeyed {
        lappend out "$u otherwise"
    }
    return [join $out ", "]
}

proc ::helpdoc::computedSentinel_ {value} {
    # A computed-sentinel value is rendered in square brackets. Keep the token
    # list in sync with vocab(computed) in helpdoc.d/vocabularies.tcl and the
    # "computed-sentinel" template in dev-tools/input_xx.xsl.
    set v [string trim $value]
    if { $v in {from_pseudopotential from_xml from_environment internal} } {
        return "\[$v\]"
    }
    return $value
}

proc ::helpdoc::unitsGlossKeyed_ {key glossed} {
    # render a keyed entry: "<i>" -> "for index <i>"; "<lo>-<hi>" ->
    # "for indices <lo>-<hi>".
    if { [string first - $key] >= 0 } {
        return "$glossed for indices $key"
    }
    return "$glossed for index $key"
}

proc ::helpdoc::txt_ref_link {content} {
    set re_ref  {(@ref)\s+(\w+([%]\w)*)}
    set re_link {(@link)\s+([.,;:]*[\w\+-]+([.,;:][\w\+-]+)*)}
    set re "($re_ref|$re_link)"
    return [regsub -all $re $content {"\3"}]
}
proc ::helpdoc::txt_tag_enter {tree node tag attr content depth} {
    variable txtDepth
    variable indentNum
    variable fid
    variable arr
    variable vargroup
    variable dimensiongroup
    variable colgroup
    variable rowgroup
    variable card
    variable mode
    variable rows
    variable cols
    variable info
    variable options
    variable options_first
    
    if { [info exists arr] } {
	unset arr
    }

    set content [formatString [trimEmpty [txt_atTags [txt_ref_link $content]]]]
    attr2array_ arr $attr

    global sourcedir
    source [file join $sourcedir txt_enter.tcl]
}


proc ::helpdoc::txt_tag_leave {tree node tag attr content depth} {
    variable fid 
    variable txtDepth   
    variable vargroup
    variable dimensiongroup
    variable colgroup
    variable rowgroup
    variable mode
    variable card
    variable rows
    variable cols
    variable arr
    variable options
    variable options_first

    attr2array_ arr $attr
    global sourcedir
    source [file join $sourcedir txt_leave.tcl]
}


proc ::helpdoc::txt_subtree {tree node newMode} {
    variable mode

    lappend mode $newMode

    set newTree [::struct::tree]
    $newTree deserialize [$tree serialize $node]

    $newTree walkproc [$newTree rootname] -order both txt_subtree_print
    $newTree destroy	

    ::tclu::lpop mode
}


proc ::helpdoc::txt_subtree_print {tree node action} {
    set depth [$tree depth $node]

    set tag        [$tree get $node tag]
    set attributes [getFromTree $tree $node attributes]
    set content    [getFromTree $tree $node text]
    
    set content [formatString [trimEmpty [txt_atTags [txt_ref_link $content]]]]

    txt_tag_${action} $tree $node $tag $attributes $content [expr $depth - 1]
}


proc ::helpdoc::printableVarDescription {tree node} {
    variable mode

    # Purpose: the description of variable in the card is printed only
    # when at least one of info, status or see records is present.

    set Info   [getTextFromDescendant $tree $node info]
    set Status [getTextFromDescendant $tree $node status]
    set See    [getTextFromDescendant $tree $node see]
    set Opt    [getTextFromDescendant $tree $node opt]

    if { ! [::tclu::lpresent $mode card] || ($Info != "" || $Status != "" || $See != "" || $Opt != "") } {
	return 1
    } 

    return 0
}
