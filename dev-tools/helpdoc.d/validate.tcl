#
# validate.tcl -- helpdoc "--strict" validation pass
#
# Walks the parsed .def tree and enforces the modernized .def conventions
# that the schema language itself cannot express. Only invoked with the
# -strict flag. Closed vocabularies are loaded from helpdoc.d/vocabularies.tcl;
# if absent the validator skips the vocabulary checks instead of crashing.
#

namespace eval ::helpdoc {
    # collected diagnostics for the current file
    variable strictErrors {}
}


proc ::helpdoc::strictErr_ {node msg} {
    # record one strict-mode validation error
    variable strictErrors
    variable tree

    set where ""
    catch {
	set tag [getFromTree $tree $node tag]
	set attr [getFromTree $tree $node attributes]
	attr2array_ _a $attr
	if { [info exists _a(name)] } {
	    set where " <$tag name=\"$_a(name)\">"
	} else {
	    set where " <$tag>"
	}
    }
    lappend strictErrors "${where}: $msg"
}


proc ::helpdoc::loadVocabularies_ {} {
    # Load the closed vocabularies into the ::helpdoc::vocab array from
    # helpdoc.d/vocabularies.tcl. Returns 1 on success, 0 if the file is
    # missing/unusable (vocabulary checks are then skipped).
    variable dir
    variable vocab

    # the vocabulary file location can be overridden via the
    # HELPDOC_VOCABULARIES environment variable (used by the test harness)
    if { [info exists ::env(HELPDOC_VOCABULARIES)] && $::env(HELPDOC_VOCABULARIES) ne {} } {
	set vfile $::env(HELPDOC_VOCABULARIES)
    } else {
	set vfile [file join $dir vocabularies.tcl]
    }

    if { ! [file exists $vfile] } {
	puts stderr "
### \[helpdoc --strict\] WARNING: vocabulary file not found: $vfile
###   The closed-vocabulary checks will be skipped.
"
	return 0
    }

    if { [catch { source $vfile } errMsg] } {
	puts stderr "
### \[helpdoc --strict\] WARNING:
###   failed to load vocabulary file:
###     $vfile
###   error: $errMsg
###   The closed-vocabulary checks will be skipped.
"
	return 0
    }

    # sanity-check the expected keys; tolerate missing ones individually
    foreach key {dim_base units_base kind computed dimensionless} {
	if { ! [info exists vocab($key)] } {
	    puts stderr "### \[helpdoc --strict\] WARNING: vocabulary '$key' missing from $vfile"
	    set vocab($key) {}
	}
    }
    return 1
}


proc ::helpdoc::parseKeyedClist_ {node text} {
    # Parse a keyed comma-list (content of <dimensionality>/<units>): each
    # token is "term", "<i>:term" or "<lo>-<hi>:term". Returns a list of
    # {lo hi term} triplets (unkeyed token has lo=hi={}).
    set result {}
    foreach tok [split $text ,] {
	set tok [string trim $tok]
	if { $tok eq {} } { continue }
	if { [string first : $tok] < 0 } {
	    # unkeyed token: default for all indices
	    lappend result [list {} {} $tok]
	    continue
	}
	set key  [string trim [lindex [split $tok :] 0]]
	set term [string trim [join [lrange [split $tok :] 1 end] :]]
	if { $term eq {} } {
	    strictErr_ $node "malformed keyed-clist token \"$tok\" (empty term)"
	    continue
	}
	if { [regexp {^([0-9]+)-([0-9]+)$} $key -> lo hi] } {
	    if { $lo > $hi } {
		strictErr_ $node "malformed keyed-clist token \"$tok\" (range lo > hi)"
		continue
	    }
	    lappend result [list $lo $hi $term]
	} elseif { [regexp {^([0-9]+)$} $key -> idx] } {
	    lappend result [list $idx $idx $term]
	} else {
	    strictErr_ $node "malformed keyed-clist token \"$tok\" (index must be <i> or <lo>-<hi>)"
	}
    }
    return $result
}




proc ::helpdoc::checkProductOfPowers_ {node expr vocabKey label {emptyToken {}}} {
    # Validate a value as a product-of-powers expression: space-separated
    # tokens, each "<base>" or "<base>^<exponent>" (exponent integer or
    # rational p/q), with each <base> in vocab($vocabKey). $label names the
    # field in diagnostics; $emptyToken, if set, is a standalone empty-product
    # literal (e.g. "dimensionless").
    variable vocab

    set haveVocab [expr {[info exists vocab($vocabKey)] && $vocab($vocabKey) ne {}}]

    foreach token [split $expr] {
        if { $token eq {} } { continue }
        # the empty-product literal is accepted as a standalone base token
        if { $emptyToken ne {} && $token eq $emptyToken } { continue }
        set base $token
        set exponent {}
        if { [string first ^ $token] >= 0 } {
            set parts [split $token ^]
            if { [llength $parts] != 2 } {
                strictErr_ $node "$label token \"$token\" is malformed (expected <base> or <base>^<exponent>)"
                continue
            }
            set base     [lindex $parts 0]
            set exponent [lindex $parts 1]
        }
        if { $base eq {} } {
            strictErr_ $node "$label token \"$token\" has an empty base"
            continue
        }
        if { $haveVocab && [lsearch -exact $vocab($vocabKey) $base] < 0 } {
            strictErr_ $node "$label base token \"$base\" is not in the closed vocabulary"
        }
        if { $exponent ne {} } {
            # exponent must be an integer, or a rational p/q with integers p,q
            if { [regexp {^[+-]?[0-9]+$} $exponent] } {
                # plain integer -- ok
            } elseif { [regexp {^([+-]?[0-9]+)/([+-]?[0-9]+)$} $exponent -> p q] } {
                if { $q == 0 } {
                    strictErr_ $node "$label token \"$token\" has a zero denominator in its exponent"
                }
            } else {
                strictErr_ $node "$label token \"$token\" has a malformed exponent \"$exponent\" (expected integer or p/q rational)"
            }
        }
    }
}


proc ::helpdoc::checkUnitsClist_ {node text} {
    # Validate a keyed comma-list of <units>; returns the triplets so the
    # caller can check index bounds.
    set triplets [parseKeyedClist_ $node $text]
    foreach t $triplets {
        checkProductOfPowers_ $node [lindex $t 2] units_base "units" dimensionless
    }
    return $triplets
}


proc ::helpdoc::checkDimClist_ {node text} {
    # Validate a keyed comma-list of <dimensionality>; returns the triplets
    # for bounds checking.
    set triplets [parseKeyedClist_ $node $text]
    foreach t $triplets {
        checkProductOfPowers_ $node [lindex $t 2] dim_base "dimensionality" dimensionless
    }
    return $triplets
}


proc ::helpdoc::checkClistBounds_ {node triplets} {
    # check that keyed-clist indices fall within the entry's start/end bounds
    variable tree

    set attr [getFromTree $tree $node attributes]
    attr2array_ a $attr

    # gather numeric bounds (var has none; dimension has start/end;
    # multidimension has clist start/end)
    set lo {}
    set hi {}
    if { [info exists a(start)] && [string is integer -strict $a(start)] } {
	set lo $a(start)
    }
    if { [info exists a(end)] && [string is integer -strict $a(end)] } {
	set hi $a(end)
    }
    # multidimension: start/end may be comma-lists -> use extremes
    if { [info exists a(start)] && [string first , $a(start)] >= 0 } {
	set vals {}
	foreach v [split $a(start) ,] {
	    set v [string trim $v]
	    if { [string is integer -strict $v] } { lappend vals $v }
	}
	if { $vals ne {} } { set lo [lindex [lsort -integer $vals] 0] }
    }
    if { [info exists a(end)] && [string first , $a(end)] >= 0 } {
	set vals {}
	foreach v [split $a(end) ,] {
	    set v [string trim $v]
	    if { [string is integer -strict $v] } { lappend vals $v }
	}
	if { $vals ne {} } { set hi [lindex [lsort -integer $vals] end] }
    }

    foreach t $triplets {
	foreach {tlo thi term} $t { break }
	if { $tlo eq {} } { continue }; # unkeyed token: no bounds to check
	if { $lo ne {} && $tlo < $lo } {
	    strictErr_ $node "keyed-clist index $tlo is below the entry's start bound ($lo)"
	}
	if { $hi ne {} && $thi > $hi } {
	    strictErr_ $node "keyed-clist index $thi is above the entry's end bound ($hi)"
	}
    }
}


proc ::helpdoc::isDimensionless_ {triplets} {
    # Returns 1 if every keyed-clist term is the literal "dimensionless",
    # 0 if any term is dimensionful, -1 if undetermined (vocabulary absent).
    variable vocab

    if { ! [info exists vocab(dimensionless)] || $vocab(dimensionless) eq {} } {
	# without the vocabulary we cannot tell; assume dimensionful so that
	# a missing <units> is flagged conservatively only when we DO know.
	return -1
    }
    if { $triplets eq {} } { return -1 }
    foreach t $triplets {
	set term [lindex $t 2]
	if { [lsearch -exact $vocab(dimensionless) $term] < 0 } {
	    return 0
	}
    }
    return 1
}


proc ::helpdoc::validateDefault_ {node} {
    # A <default> has EITHER a text body XOR <case> children; the "kind"
    # attribute must match the body shape; a conditional default's last
    # <case> must be a testless fallback.
    variable tree
    variable vocab

    set attr [getFromTree $tree $node attributes]
    attr2array_ a $attr
    set kind literal
    if { [info exists a(kind)] && [string trim $a(kind)] ne {} } {
	set kind [string trim $a(kind)]
    }

    # kind must be in vocabulary (when the vocabulary is available)
    if { [info exists vocab(kind)] && $vocab(kind) ne {} } {
	if { [lsearch -exact $vocab(kind) $kind] < 0 } {
	    strictErr_ $node "default kind \"$kind\" is not in the closed vocabulary"
	}
    }

    # gather case children and the text body
    set caseNodes {}
    foreach child [$tree children $node] {
	if { [getFromTree $tree $child tag] eq "case" } {
	    lappend caseNodes $child
	}
    }
    set body [string trim [getFromTree $tree $node text]]

    set hasBody  [expr {$body ne {}}]
    set hasCases [expr {[llength $caseNodes] > 0}]

    # EITHER text body XOR case children, never both / never neither
    if { $hasBody && $hasCases } {
	strictErr_ $node "a <default> must have EITHER a text body XOR <case> children, not both"
    }
    if { ! $hasBody && ! $hasCases } {
	strictErr_ $node "a <default> is empty: it must have either a text body or <case> children"
    }

    # kind must match the body shape
    if { $kind eq "conditional" } {
	if { ! $hasCases } {
	    strictErr_ $node "default kind=\"conditional\" requires <case> children"
	}
    } else {
	if { $hasCases } {
	    strictErr_ $node "default with <case> children must use kind=\"conditional\" (got kind=\"$kind\")"
	}
    }

    # a conditional default's last <case> must be a testless fallback
    if { $hasCases } {
	set last [lindex $caseNodes end]
	set lattr [getFromTree $tree $last attributes]
	attr2array_ la $lattr
	if { [info exists la(test)] && [string trim $la(test)] ne {} } {
	    strictErr_ $node "the last <case> of a conditional <default> must be a fallback (without a \"test\")"
	}
	# all non-last cases SHOULD carry a test
	foreach c [lrange $caseNodes 0 end-1] {
	    set cattr [getFromTree $tree $c attributes]
	    attr2array_ ca $cattr
	    if { ! [info exists ca(test)] || [string trim $ca(test)] eq {} } {
		strictErr_ $node "every <case> before the fallback must carry a \"test\" attribute"
	    }
	}
    }

    # a kind=computed default's body must be a computed-sentinel token
    if { $kind eq "computed" && $hasBody } {
	if { [info exists vocab(computed)] && $vocab(computed) ne {} } {
	    if { [lsearch -exact $vocab(computed) $body] < 0 } {
		strictErr_ $node "computed-default sentinel \"$body\" is not in the closed vocabulary"
	    }
	}
    }
}


proc ::helpdoc::optValSet_ {optionsNode} {
    # Collect the flat list of all values declared by an <options> block;
    # each <opt> "val" (and optional "alias") attribute is itself a comma
    # list. Quotes are stripped.
    variable tree

    set vals {}
    foreach child [$tree children $optionsNode] {
	if { [getFromTree $tree $child tag] ne "opt" } { continue }
	set cattr [getFromTree $tree $child attributes]
	attr2array_ ca $cattr
	foreach key {val alias} {
	    if { ! [info exists ca($key)] } { continue }
	    foreach v [split $ca($key) ,] {
		set v [string trim $v]
		set v [string trim $v {'\"}]
		set v [string trim $v]
		if { $v ne {} } { lappend vals $v }
	    }
	}
    }
    return $vals
}


proc ::helpdoc::validateEnumDefault_ {node} {
    # If an entry has both an <options> block and a literal <default>, the
    # default value must be one of the declared opt values. Non-literal
    # defaults (kind=conditional/computed) are exempt.
    variable tree

    set optionsNode {}
    set defaultNode  {}
    foreach child [$tree children $node] {
	# if-chain, not switch: "default" is a real tag name and would clash
	# with switch's wildcard fallback arm.
	set ctag [getFromTree $tree $child tag]
	if { $ctag eq "options" } {
	    set optionsNode $child
	} elseif { $ctag eq "default" } {
	    set defaultNode $child
	}
    }
    if { $optionsNode eq {} || $defaultNode eq {} } { return }

    # only a *literal* default has a single value to check
    set dattr [getFromTree $tree $defaultNode attributes]
    attr2array_ da $dattr
    set dkind literal
    if { [info exists da(kind)] && [string trim $da(kind)] ne {} } {
	set dkind [string trim $da(kind)]
    }
    if { $dkind ne "literal" } { return }

    set defText [string trim [getFromTree $tree $defaultNode text]]
    if { $defText eq {} } { return }
    set defValue [string trim [string trim $defText {'\"}]]
    # a blank literal default is a legitimate "unset" sentinel; exempt it
    if { $defValue eq {} } { return }

    set allowed [optValSet_ $optionsNode]
    if { $allowed eq {} } { return }

    if { [lsearch -exact $allowed $defValue] < 0 } {
	strictErr_ $node "enumerated <default> value \"$defValue\" is not among the declared <opt> -val/-alias choices: [join $allowed {, }]"
    }
}


proc ::helpdoc::caseNodesOf_ {node} {
    # Return the <case> child nodes of $node, in document order.
    variable tree
    set caseNodes {}
    foreach child [$tree children $node] {
	if { [getFromTree $tree $child tag] eq "case" } {
	    lappend caseNodes $child
	}
    }
    return $caseNodes
}


proc ::helpdoc::validateConditionalCases_ {node label caseNodes} {
    # The last <case> must be a testless fallback; every preceding <case>
    # must carry a "test" attribute. $label names the element in diagnostics.
    variable tree

    if { [llength $caseNodes] == 0 } { return }
    set last [lindex $caseNodes end]
    set lattr [getFromTree $tree $last attributes]
    attr2array_ la $lattr
    if { [info exists la(test)] && [string trim $la(test)] ne {} } {
	strictErr_ $node "the last <case> of a conditional <$label> must be a fallback (without a \"test\")"
    }
    foreach c [lrange $caseNodes 0 end-1] {
	set cattr [getFromTree $tree $c attributes]
	attr2array_ ca $cattr
	if { ! [info exists ca(test)] || [string trim $ca(test)] eq {} } {
	    strictErr_ $node "every <case> before the fallback of a conditional <$label> must carry a \"test\" attribute"
	}
    }
}


proc ::helpdoc::validateConditionalUnitsDim_ {node label} {
    # Validate a <dimensionality> or <units> element ($label), in both the
    # plain text-body form (content checked by the caller) and the
    # conditional form (kind="conditional" with <case> children). Returns
    # {isConditional allDimensionless}: allDimensionless is 1/0/-1 and
    # meaningful only for dimensionality.
    variable tree
    variable vocab

    set attr [getFromTree $tree $node attributes]
    attr2array_ a $attr
    set kind plain
    if { [info exists a(kind)] && [string trim $a(kind)] ne {} } {
	set kind [string trim $a(kind)]
    }

    set caseNodes [caseNodesOf_ $node]
    set body [string trim [getFromTree $tree $node text]]
    set hasBody  [expr {$body ne {}}]
    set hasCases [expr {[llength $caseNodes] > 0}]

    if { ! $hasCases } {
	# plain text-body form: caller handles content validation
	if { $kind eq "conditional" } {
	    strictErr_ $node "<$label> kind=\"conditional\" requires <case> children"
	}
	return {0 -1}
    }

    # conditional form
    if { $hasBody } {
	strictErr_ $node "a conditional <$label> must have <case> children XOR a text body, not both"
    }
    if { $kind ne "conditional" } {
	strictErr_ $node "<$label> with <case> children must use kind=\"conditional\" (got kind=\"$kind\")"
    }

    validateConditionalCases_ $node $label $caseNodes

    # validate each case's value text as a product-of-powers expression
    set vocabKey [expr {$label eq "units" ? "units_base" : "dim_base"}]
    set emptyTok dimensionless
    set allDimless 1
    set anyKnown 0
    foreach c $caseNodes {
	set ctext [string trim [getFromTree $tree $c text]]
	if { $ctext eq {} } {
	    strictErr_ $node "a <case> of a conditional <$label> has an empty value"
	    continue
	}
	checkProductOfPowers_ $c $ctext $vocabKey $label $emptyTok
	if { $label eq "dimensionality" } {
	    if { [info exists vocab(dimensionless)] && $vocab(dimensionless) ne {} } {
		set anyKnown 1
		if { [lsearch -exact $vocab(dimensionless) $ctext] < 0 } {
		    set allDimless 0
		}
	    }
	}
    }

    if { $label eq "dimensionality" } {
	if { ! $anyKnown } { set allDimless -1 }
	return [list 1 $allDimless]
    }
    return {1 -1}
}


proc ::helpdoc::entryType_ {node} {
    # Return the upper-cased data type of an entry, from its own "-type"
    # attribute or inherited from an enclosing typed *group. Returns "" when
    # undetermined (no type, or an ambiguous comma-list group type).
    variable tree

    set attr [getFromTree $tree $node attributes]
    attr2array_ a $attr
    if { [info exists a(type)] && [string trim $a(type)] ne {} } {
	set t [string trim $a(type)]
	if { [string first , $t] >= 0 } { return {} }
	return [string toupper $t]
    }

    # no own -type: look one level up for an enclosing typed group
    set parent [$tree parent $node]
    if { $parent ne {} } {
	switch -glob -- [getFromTree $tree $parent tag] {
	    *group {
		set pattr [getFromTree $tree $parent attributes]
		attr2array_ pa $pattr
		if { [info exists pa(type)] && [string trim $pa(type)] ne {} } {
		    set t [string trim $pa(type)]
		    if { [string first , $t] >= 0 } { return {} }
		    return [string toupper $t]
		}
	    }
	}
    }
    return {}
}


proc ::helpdoc::validateVarEntry_ {node} {
    # enforce rules 1, 2, 4 for a single var/dimension-class entry
    variable tree

    set tag [getFromTree $tree $node tag]

    # A childless entry is a bare re-reference of an already-defined
    # variable; it has no body and cannot declare <dimensionality>/<units>.
    if { [llength [$tree children $node]] == 0 } {
	return
    }

    # collect the relevant immediate children
    set dimNode   {}
    set unitsNode {}
    foreach child [$tree children $node] {
	switch -exact -- [getFromTree $tree $child tag] {
	    dimensionality { set dimNode   $child }
	    units          { set unitsNode $child }
	}
    }

    # <dimensionality>/<units> may be inherited from an enclosing *group
    if { $dimNode eq {} || $unitsNode eq {} } {
	set parent [$tree parent $node]
	if { $parent ne {} } {
	    switch -glob -- [getFromTree $tree $parent tag] {
		*group {
		    foreach child [$tree children $parent] {
			switch -exact -- [getFromTree $tree $child tag] {
			    dimensionality {
				if { $dimNode eq {} } { set dimNode $child }
			    }
			    units {
				if { $unitsNode eq {} } { set unitsNode $child }
			    }
			}
		    }
		}
	    }
	}
    }

    # rule 1: <dimensionality> is REAL-only -- required for a REAL entry,
    # forbidden for a non-REAL one, skipped when the type is undetermined.
    set etype [entryType_ $node]
    if { $etype eq "REAL" } {
	if { $dimNode eq {} } {
	    strictErr_ $node "$tag entry is REAL but has no <dimensionality> (mandatory under --strict)"
	    return
	}
    } elseif { $etype ne {} } {
	# a determined non-REAL type
	if { $dimNode ne {} } {
	    strictErr_ $node "$tag entry has -type $etype but carries a <dimensionality> (forbidden under --strict: <dimensionality> is REAL-only)"
	}
	return
    } else {
	# type cannot be determined: neither require nor forbid
	if { $dimNode eq {} } { return }
    }

    # validate the <dimensionality> content (plain or conditional form)
    foreach {dimCond dimAllDimless} [validateConditionalUnitsDim_ $dimNode dimensionality] { break }
    if { $dimCond } {
	set dimless $dimAllDimless
    } else {
	set dimText [string trim [getFromTree $tree $dimNode text]]
	set dimTriplets [checkDimClist_ $dimNode $dimText]
	checkClistBounds_ $node $dimTriplets
	set dimless [isDimensionless_ $dimTriplets]
    }

    # rule 2: units mandatory whenever dimensionality is dimensionful
    if { $dimless == 0 } {
	if { $unitsNode eq {} } {
	    strictErr_ $node "$tag entry is dimensionful but has no <units>"
	}
    }

    # validate the <units> content (if present), plain or conditional form
    if { $unitsNode ne {} } {
	foreach {unitsCond unitsDummy} [validateConditionalUnitsDim_ $unitsNode units] { break }
	if { ! $unitsCond } {
	    set unitsText [string trim [getFromTree $tree $unitsNode text]]
	    set unitsTriplets [checkUnitsClist_ $unitsNode $unitsText]
	    checkClistBounds_ $node $unitsTriplets
	}
    }
}


proc ::helpdoc::validateRequiredNoDefault_ {node} {
    # An entry whose <status> reads "required" must not also carry a
    # <default> child.
    variable tree

    set statusRequired 0
    set defaultNode {}
    foreach child [$tree children $node] {
	set ctag [getFromTree $tree $child tag]
	if { $ctag eq "status" } {
	    set stext [string trim [getFromTree $tree $child text]]
	    set stext [string trim $stext "{} \t\n"]
	    if { [string equal -nocase $stext "required"] } {
		set statusRequired 1
	    }
	} elseif { $ctag eq "default" } {
	    set defaultNode $child
	}
    }

    if { $statusRequired && $defaultNode ne {} } {
	strictErr_ $node "entry has <status> required but also carries a <default> (forbidden under --strict: a required entry has no default)"
    }
}


proc ::helpdoc::strictValidate {file} {
    # Top-level --strict pass: walk the parsed tree of $file enforcing the
    # modernized .def conventions; abort if any violation is found.
    variable tree
    variable strictErrors

    puts "\n***\n*** \[helpdoc --strict\] validating $file\n***\n"

    set strictErrors {}
    loadVocabularies_

    foreach node [$tree descendants root] {
	set tag [getFromTree $tree $node tag]
	# if-chain, not switch: "default" is a real tag name (see above).
	if { $tag eq "var" || $tag eq "dimension" || $tag eq "multidimension" \
	     || $tag eq "col" || $tag eq "row" } {
	    validateVarEntry_ $node
	    validateEnumDefault_ $node
	    validateRequiredNoDefault_ $node
	} elseif { $tag eq "default" } {
	    validateDefault_ $node
	}
    }

    if { [llength $strictErrors] > 0 } {
	puts stderr "
### \[helpdoc --strict\] VALIDATION FAILED for $file
### [llength $strictErrors] violation(s):
"
	foreach e $strictErrors {
	    puts stderr "  *$e"
	}
	puts stderr ""
	::tclu::abort "strict validation failed for $file"
    }

    puts "*** \[helpdoc --strict\] OK - $file passed strict validation\n"
}
