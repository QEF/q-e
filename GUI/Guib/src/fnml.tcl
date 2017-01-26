#
# $RCSfile: fnml.tcl,v $ --
#
#      This file contains a Tcl parser for the Fortran
#      namelist. Parsed namelist is returned as "variable value"
#      pairs suitable for the "array set" command. Implicit assigments
#      of one-dimensional arrays are returned as explicit assigments,
#      except the "r*c" and "r*" forms are left intact. Implicit
#      assignment for multi-dimensional arrays is not supported and
#      results in parsing error.
#    
#
# Copyright (c) 2017  Anton Kokalj   Email: tone.kokalj@ijs.si
#
#
# This file is distributed under the terms of the GNU General Public
# License. See the file `COPYING' in the root directory of the present
# distribution, or http://www.gnu.org/copyleft/gpl.txt .
#
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# ANTON KOKALJ BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
# AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
#
# $Id: tclUtils.tcl,v 1.1.1.1 2014/09/08 12:15:21 tone Exp $ 
#


namespace eval ::fnml {
    variable errInfo {}
    variable errCode PARSE_ERROR
}


#****f* ::fnml/::fnml::parseContent
# SYNOPSIS
proc ::fnml::parseContent {nmlName content} {    
    # PURPOSE    
    # Parse the content of the namelist and return the parsed namelist
    # as variable-value pairs in format suitable for "array get"
    # command. Note that implicit array assigments are returned as
    # explicit assigments, except the "r*c" and "r*" forms are left intact.
    #    
    # ARGUMENTS
    # * nmlName -- name of the namelist
    # * content -- content of the namelist w/o the namelist's begin and end tags (i.e. &name, / or &end)
    #
    # RETURN VALUE
    # Returns the parsed namelist as variable-value pairs in format
    # suitable for "array get" command.
    #
    # ERRORS
    # Upon detecting parsing (syntax) error in the namelist an error
    # is triggered and error-code is set to PARSE_ERROR
    #
    # SOURCE
    variable errInfo 
    variable errCode 

    set nmlName &[string trimleft $nmlName &]
    
    set err_prefix "error parsing \"$nmlName\" namelist"
    
    # split namelist-content by delimiters [\t\n, ] and remove comments
    
    set record [split [::fnml::preparse_ $content] =]

    if { [catch {set nr [llength $record]}] } {
	# if [llength] fails there is an "unmatched open quote in list" --> syntax error in namelist
	error "$err_prefix: unmatched open quote in the namelist" $errInfo $errCode
    }
    if { $nr == 1 } {
	# if nr == 1 --> there is no "=" in the namelist, but the namelist isn't empty
	error "$err_prefix: no = operator found in a non-empty namelist" $errInfo $errCode
    }	

    foreach field $record {

	incr ind	    
	
	if { ! [info exists next_var] } {

	    # namelist's first field
	    
	    if { $ind == 1 && ([catch {set nf_ [llength $field]}] || $nf_ > 1) } {
		# the first field in the namelist should be a single word,
		# but here we have something like {var1 var2 ... = ...}
		error "$err_prefix: problem in the specs of the first variable in the namelist" $errInfo $errCode
	    }	    
	    if { [catch {set next_var [lindex $field end]}] } {
		error "$err_prefix: unmatched open quote in the namelist" $errInfo $errCode
	    }
	    if { $next_var eq {} } {
		error "$err_prefix: namelist data probably start with the = operator" $errInfo $errCode
	    } 
	} else {
	    set var $next_var

	    if { $var == {} } {
		error "$err_prefix: probably two consecutive = operators" $errInfo $errCode
	    }
	    
	    if { $nr == $ind } {
		# last field
		set values $field
	    } else {
		if { [catch {set nf [llength $field]}] } {
		    error "$err_prefix: probably unmatched open quote in the namelist" $errInfo $errCode
		}
		set values   [lrange $field  0  $nf-2]
		set next_var [lindex $field end]
	    }

	    #
	    # assing value(s) to variable
	    #
	    
	    if { [llength $values] <= 1 } {
		#
		# scalar variable (or explicitly specified element of an array)
		#
		
		set nml($var) [string trim $values "\{\}\t\n "]

	    } else {
		#
		# it's an array
		#
		set index 1
		
		# fully implicit spec, i.e "var" instead of "var(index)"
		set re {\w+[%\w]*}
		
		# regexp for "post-index" (:n) form (note the index "n" of this form is ignored by gfortran and ifort)
		set re_post {(\w+[%\w]*)\( *: *([+-]?[0-9]*) *\)}
    
		# regexp for "pre-index" (n), (n:), or (n:m) forms, which all implies "n"
		set re_pre {(\w+[%\w]*)\( *([+-]?[0-9]+) *:? *([+-]?[0-9]*) *\)}
		
		# explicit multidimensional index specs, i.e. (m,n,k)
		set rem_expl {(\w+[%\w]*)\( *((([+-]?[0-9]+) *,)+( *[+-]?[0-9]+)) *\)}

		# implicit multidimensional index specs, various forms of (n:m,k:l*): this is NOT SUPPORTED
		set rem_impl {(\w+[%\w]*)\( *([+-]?[0-9]*) *: *([+-]?[0-9]*) *,}
		
		if { [regexp $rem_expl $var] || [regexp $rem_impl $var] } {
		    
		    error "$err_prefix: error while parsing variable $var of \"$nmlName\" namelist: implicit indexing for assignment of multidimensional arrays is not supported. Use explicit assigment instead." $errInfo $errCode
		    
		} elseif { [regexp $re_pre $var] } {
		    
		    set index [regsub $re_pre $var {\2}]
		    set var   [lindex [split $var \(] 0]		    

		} elseif { [regexp $re_post $var] || [regexp $re $var] } {
		    
		    set index 1
		    set var   [lindex [split $var \(] 0]
		} else {
		    error "shouldn't happen: something wrong when parsing the namelist $nmlName. Please report the bug along with the namelist" $errInfo $errCode
		}

		foreach val $values {
		    set name ${var}($index)
		    set nml($name) [string trim $val "\{\}\t\n "]
		    incr index		
		}
	    }
	}
    }

    return [array get nml]
}
#******


#****f* ::fnml/::fnml::parseContent
# SYNOPSIS
proc ::fnml::getContent {nmlName fileChannel {match_first 1} {delim_re {[, \t\n]}}} {
    # PURPOSE
    # Extract the content of the specified namelist from the
    # file-channel. The content of the namelist is trimmed from the
    # namelist's begin and end tags (i.e. &name, / or &end).
    #    
    # ARGUMENTS
    # * nmlName -- name of the namelist
    # * fileChannel -- file-channel
    # * match_first -- [optional] if true then the first non-empty line must start with the $nmlName word (default = 1)
    # * delim_re -- [optional] regexp for the allowed delimeter characters (default = {[, \t\n]})
    #
    # RETURN VALUE
    # The content of the namelist w/o the namelist's begin and end tags (i.e. &name, / or &end).
    #
    # ERRORS
    # Upon a parsing error (either EOF or if the namelist was not
    # found) an error is triggered and error-code is set to PARSE_ERROR
    #
    # SOURCE
    variable errInfo 
    variable errCode 

    set nmlName &[string trimleft $nmlName &]

    while {1} {
	incr nl
	set line [readline_ $fileChannel]
	scan $line %s name
	if { $nl == 1 } { set name_1st $name }
	
	if { [string toupper $name] eq [string toupper $nmlName] } {
	    if { $match_first && $nl > 1 } {
		# the namelist wasn't found at the first non-empty line 
		error "expected namelist $nmlName, but got \"$name_1st\""  $errInfo $errCode
	    }
	    break
	}
	if { [eof $fileChannel] } {
	    error "end of file while searching for namelist $nmlName"  $errInfo $errCode
	}
    }    
    
    # strip the nmlName from line		
    set line [regsub -nocase "^$nmlName" $line {}]
        
    set s ';  # single-quote
    set d \"; # double-quote    
    set q -1;       # stores the quoted character
    set inside_q 0; # inside|outside the quoted string
    set result  {}; # resulting correctly-delimiter-split list
    set end      0; # will be one when encountering end-of-namelist

    while {1} {

	foreach c [split $line\n {}] {
	    
	    if { ! $inside_q } {
		
		# outside quoted string

		if { $c eq "!" } {
		    # it's comment, skip the rest of the line
		    append result \n
		    break
		}
		if { $c eq "/" } {
		    # it's the end of namelist character
		    set end 1
		    break
		}

		if { $c eq $s || $c eq $d } {
		    # string starting quote
		    set q $c
		    set inside_q 1
		}
		    
		if { [regexp $delim_re $c] } {
		    if { [info exists word]} {
			# delimiter outside the string --> it's delimiter
			if { $word eq "&end" } {
			    #  it's the end of namelist character
			    set end 1
			    break
			}
			unset word
		    }
		} else {
		    append word $c
		}		
	    } else {
		# inside quoted string
	    
		if { $c eq $q } {
		    # string ending quote
		    set inside_q 0
		}
		append word $c
	    }
	    append result $c
	}
	
	if { $end } {
	    return $result
	}
	if { [eof $fileChannel] && $end == 0 } {
	    error "end of file, while reading namelist $nmlName"  $errInfo $errCode
	}

	set line [readline_ $fileChannel]
    }
}
#******



#------------------------------------------------------------------------
# read a next non-empty line from a file channel and return the line trimed by "spaces"
# --
proc ::fnml::readline_ {fileChannel} {
    while {[gets $fileChannel line] > -1 && [regexp {^\s*$} $line] == 1 } {}
    return [string trim $line "\n\t "]
}


#------------------------------------------------------------------------
# This proc pre-parse the namelist content and performs 3 operations:
#
# 1. removes comments
#
# 2. read just up to the end-of-namelist and trims the rest
#
# 3. split the namelist content by delimiter ignoring the delimiter
#    inside the fortran strings and return the content as list. The
#    effect of this operation is that any fortran string is returned
#    as a single list element and within the string the delimiter,
#    end-of-namelist, and comment-character have no effect, e.g. the
#    latter does not split the string into two (or more) list
#    elements.
#
# Default value of delimeter is {[, \t\n]
#--
proc ::fnml::preparse_ {content {delim_re {[, \t\n]}}} {

    set content [string trim $content "\n\t "]
    
    set s ';  # single-quote
    set d \"; # double-quote
    
    set q -1;       # stores the quoted character
    set inside_q 0; # inside|outside the quoted string
    set inside_p 0; # inside|outside the parenthesis (i.e. for array specs and complex specs)
    set result {};  # resulting correctly-delimiter-split list
    set comment  0; # is there a comment 0|1
    
    foreach c [split $content {}] {
	if { !$inside_q } {
	    # outside quoted string
	    if { $c eq "\n" } {
		set comment 0
	    }
	    if { $c eq "!" } {
		set comment 1
		continue
	    }
	    if { $c eq "/" } {
		# it's the end of namelist character
		break
	    }
	    if { $c eq "(" } {
		set inside_p 1
	    }
	    if { $c eq ")" } {
		set inside_p 0
	    }
	    if { ! $comment } {
		if { $c eq $s || $c eq $d } {
		    # string starting quote
		    set q $c
		    set inside_q 1
		}
	    
		if { [regexp $delim_re $c] || $c eq "=" } {
		    # N,B,: the "=" is also a delimiter, but of special kind, i.e., it is used to
		    # delimit variableName from value, hence it needs a special treatment in var=value
		    # case (i.e. if "=" is no surronded by spaces)
		    
		    if { [info exists word]} {
			# delimiter outside the string --> it's delimiter
			if { $word eq "&end" } {
			    #  it's the end of namelist character
			    return $result
			}
			if { ! $inside_p } {
			    lappend result $word
			    if { $c eq "=" } { lappend result = }
			    unset word
			} else {
			    append word $c
			}
		    } elseif { $c eq "=" } {
			# this can probably happen in case of syntax error ... if so fnml::parseContent will report an error
			lappend result =
		    }
		} else {
		    append word $c
		}
	    }
	} else {
	    # inside quoted string
	    
	    if { $c eq $q } {
		# string ending quote
		set inside_q 0
	    }
	    append word $c
	}
    }
    if { [info exists word] && $word ne {} && $word ne "&end" } {
	lappend result $word
    }
    #puts "preparse_:  $result"
    return $result
}
