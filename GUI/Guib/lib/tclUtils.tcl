#
# $RCSfile: tclUtils.tcl,v $ --
#
#      This file contains the Tone Kokalj's Tcl utilities functions.
#
# Copyright (c) 2003  Anton Kokalj   Email: tone.kokalj@ijs.si
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
# $Id: tclUtils.tcl,v 1.4 2004-03-17 10:09:33 kokalj Exp $ 
#

#------------------------------------------------------------------------
#****h* TclLib/::tclu
#  NAME
#    TCLU == Tone Kokalj's Tcl Utilities
#                          ^^^ ^ == tclu
#  COPYRIGHT
#    2001--2004 (c) by Tone Kokalj
#  AUTHOR
#    Tone Kokalj
#  CREATION DATE
#    Starting on Fri Dec 14 08:32:59 CET 2001
#  MODIFICATION HISTORY
#    Writing a few functions (based on XCRYSDEN's Tcl functions)
#  NOTES
#    Here are the basic Tcl-only utility functions
#******
#------------------------------------------------------------------------

package provide tclu 0.1

namespace eval ::tclu:: {
    variable DEBUG      0
    variable DEBUG_FILE 0
    variable debug
    variable error
    variable public  
    variable private 
    variable ifexists
    variable _line "------------------------------------------------------------------------"
    variable tempFile
    variable nonblocking
    
    set tempFile(counter) 0
    set tempFile(files)   {}
    set tempFile(dirs)    {}
    set nonblocking(counter) -1

    namespace export usage
    namespace export writeFile
    namespace export readFile
    namespace export DEBUG
    namespace export ERROR
    namespace export errorDialog
    namespace export putsFlush
    namespace export printf
    namespace export printfFlush
    namespace export fprintf
    namespace export fprintfFlush
    namespace export public
    namespace export getPublic
    namespace export private
    namespace export getPrivate
    namespace export lsearchCase
    namespace export range
    namespace export newset
    namespace export conset
    namespace export ifincr
    namespace export scan
    namespace export removeElementFromList
    namespace export addElementToList
    namespace export lindexEnd
    #namespace export format
    namespace export tempFile
    namespace export tempDir
    namespace export stringMatch
    namespace export nonblocking
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::usage
#  NAME
#    ::tclu::usage -- prints a USAGE message 
#  USAGE
#    ::tclu::usage ?channelID? text
#
#  DESCRIPTION
#    It is used for printing the usage of procs
#  ARGUMENTS
#    args  - comprises a ?channelID? and text, where channelID is a 
#            channel identifier such as returned from call to open or 
#            socket, and text is a usage-text to be printed
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::usage stderr "::tclu::usage ?channelID? text"
#********
#------------------------------------------------------------------------

proc ::tclu::usage {args} {
    variable _line

    if { [llength $args] > 2 } {
	::tclu::usage "::tclu::usage ?channelID? text"
    }
    set channel stderr
    set text    [lindex $args end]
    if { [llength $args] == 2 } {
	set channel [lindex $args 0]
    }
    puts $channel [format $_line\nUsage: $text\n$_line\n]
    flush $channel
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::writeFile
#  NAME
#    ::tclu::writeFile -- write a content to a file
#  USAGE
#    ::tclu::writeFile filename content ?flag?
#
#  DESCRIPTION
#    Writes a content to a file.
#  ARGUMENTS
#    filename   - name of the file to be written
#    content    - content to be writtento a file
#    flag       - access flag to "open" command (i.e. w, w+, a, a+)
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::writeFile test.txt {This is a test.txt file} w
#********
#------------------------------------------------------------------------

proc ::tclu::writeFile {filename content {flag w}} {
    set   fID  [open $filename $flag]
    puts  $fID $content
    flush $fID
    close $fID
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::readFile
#  NAME
#    ::tclu::readFile -- read a file and returns its content
#  USAGE
#    ::tclu::readFile ?-nonewline? filename
#
#  DESCRIPTION
#    Reads a file and returns its content.
#  ARGUMENTS
#    args   - comprises a ?-nonewline? option and filename
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::readFile -nonewline test.dat
#********
#------------------------------------------------------------------------

proc ::tclu::readFile {args} {
    set filename [lindex $args end]
    set fID      [open $filename r] 
    set bytes    [file size $filename]
    if { [llength $args] == 2 } {
	set flag   [lindex $args 0]
	set output [read $flag $fID]
    } else {
	set output [eval read $fID $bytes]
    }
    close  $fID
    return $output
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::DEBUG
#  NAME
#    ::tclu::DEBUG -- prints a DEBUG message
#  USAGE
#    ::tclu::DEBUG ?-stderr? message
#
#  DESCRIPTION
#    Prints a DEBUG message.
#  ARGUMENTS
#    args   - comprises a ?-stderr? and message
#  RETURN VALUE
#    Prints a debug message.
#  EXAMPLE
#    ::tclu::DEBUG -stderr "Variable A has a value: $A"
#  TODO
#    It would be good if debug would also print a message to a DEBUG file,
#    because some programs writes its output to stdout and the DEBUG 
#    messages would interfere that. So the DEBUG file would circumvent that
#********
#------------------------------------------------------------------------

proc ::tclu::DEBUG {args} {
    variable DEBUG 
    variable DEBUG_FILE
    variable debug

    if { [lindex $args 0] == "-stderr" } {
	set channel stderr
	set text [lrange $args 1 end]	
    } else {
	set channel stdout
	set text $args
    }

    if { $DEBUG } {	
	puts $channel "DEBUG::   $text"
	flush $channel
    }

    if { $DEBUG_FILE } {
	if { ! [info exists debug(channel)] } {
	    set debug(channel) [open DEBUG w]
	    set debug(id) 0
	}
	if { [info level] > 1 } {
	    puts $debug(channel) "--\nDEBUG message \#$debug(id); Call from: [info level -1]\n$text"
	} else {
	    puts $debug(channel) "--\nDEBUG message \#$debug(id); Call from top-level\n$text"
	}   
	flush $debug(channel)
	incr debug(id)
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::ERROR
#  NAME
#    ::tclu::ERROR -- generates an error message
#  USAGE
#    ::tclu::ERROR errMsg ?info? ?code?
#
#  DESCRIPTION
#    Generates an ERROR message.
#  ARGUMENTS
#    errMsg  - error message
#    info    - used to initialize the global variable errorInfo, 
#              which is used to accumulate a stack trace
#    error   - its value is stored in Tcl's errorCode variable
#  RETURN VALUE
#    Generates and error message.
#  EXAMPLE
#    ::tclu::ERROR {an error hac occurred because variable xyz is not define}
#********
#------------------------------------------------------------------------

proc ::tclu::ERROR {errMsg {info {}} {code {}}} {    
    variable error
    variable _line
    set error(errMsg) "\n$_line\nERROR: $errMsg\n$_line\n"
    set error(info) $info
    set error(code) $code
    uplevel 1 {
	global ::tclu::error
	error $::tclu::error(errMsg) $::tclu::error(info) $::tclu::error(code)
    }
}



#------------------------------------------------------------------------
#****f* ::tclu/::tclu::errorDialog
#  NAME
#    ::tclu::errorDialog -- prints an error message
#  USAGE
#    ::tclu::errorDialog errMsg
#
#  DESCRIPTION
#    Prints an errorDialog message either to stderr or makes a call to
# ::tku::errorDialog when Tk package is present
#
#  ARGUMENTS
#    errMsg  - error message
#  RETURN VALUE
#    Returns "-code returns", hence the caller proc will return.
#  EXAMPLE
#    ::tclu::errorDialog {an error hac occurred because variable xyz is not define}
#********
#------------------------------------------------------------------------

proc ::tclu::errorDialog {errMsg} {
    if { [catch {package present Tk}] } {
	puts stderr "ERROR: $errMsg"
	flush stderr
    } else {
	tk_messageBox -title ERROR -message "ERROR: $errMsg" -type ok -icon error
    }
    return -code return ""
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::putsFlush
#  NAME
#    ::tclu::putsFlush -- Tcl "puts" + "flush"
#  USAGE
#    ::tclu::putsFlush ?-nonewline? ?channelId? string
#
#  DESCRIPTION
#    Identical to Tcl's puts, but invoke the flush immediately after.
#    See puts man-page of Tcl.
#********
#------------------------------------------------------------------------
			
proc ::tclu::putsFlush {args} {
    # puts ?-nonewline? ?channelId? string
    set ind 0
    set flags "" 
    if { [lindex $args $ind] == "-nonewline" } {
	set flags "-nonewline"
	incr ind
    }
    if { [llength [lrange $args $ind end]] == 1 } {
	set channel stdout
    } else {
	set channel [lindex $args $ind]
	incr ind
    }
    eval puts $flags $channel [lrange $args $ind end]
    flush $channel
}

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::printf
#  NAME
#    ::tclu::printf -- a C-style printf
#  USAGE
#    ::tclu::printf format ?arg? ?arg? ...
#
#  DESCRIPTION
#    Prints to stdout in like C-style printf.
#  ARGUMENTS
#    format  - format description
#    args    - arguments to print
#  RETURN VALUE
#    Prints arguments to stdout in C-style printf.
#  EXAMPLE
#    ::tclu::printf "%s%d\n" "Value of i:" $i
#********
#------------------------------------------------------------------------

proc ::tclu::printf {format args} {
    puts stdout [eval format [list $format] $args]
}

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::printfFlush
#  NAME
#    ::tclu::printfFlush -- a C-style printf with stdout flush
#  USAGE
#    ::tclu::printfFlush format ?arg? ?arg? ...
#
#  DESCRIPTION
#    Prints to stdout in like C-style printf and flushes stdout.
#  ARGUMENTS
#    format  - format description
#    args    - arguments to print
#  RETURN VALUE
#    Prints arguments to stdout in C-style printf and fulshes stdout.
#  EXAMPLE
#    ::tclu::printfFlush "%s%d\n" "Value of i:" $i
#********
#------------------------------------------------------------------------
			
proc ::tclu::printfFlush {format args} {
    puts stdout [eval format [list $format] $args]
    flush stdout
}

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::fprintf
#  NAME
#    ::tclu::fprintf -- a C-style fprintf
#  USAGE
#    ::tclu::fprintf channel format ?arg? ?arg? ...
#
#  DESCRIPTION
#    Prints to $channel in like C-style fprintf.
#  ARGUMENTS
#    channel - a channel identifier such as returned from call to 
#              open or socket
#    format  - format description
#    args    - arguments to print
#  RETURN VALUE
#    Prints arguments to $channel in C-style fprintf.
#  EXAMPLE
#    ::tclu::printf $channelID "%s%d\n" "Value of i:" $i
#********
#------------------------------------------------------------------------

proc ::tclu::fprintf {channel format args} {
    puts $channel [eval format [list $format] $args]
}

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::fprintfFlush
#  NAME
#    ::tclu::fprintf -- a C-style fprintf + flush of channel
#  USAGE
#    ::tclu::fprintf channel format ?arg? ?arg? ...
#
#  DESCRIPTION
#    Prints to $channel in like C-style fprintf and flashes the $channel
#  ARGUMENTS
#    channel - a channel identifier such as returned from call to 
#              open or socket
#    format  - format description
#    args    - arguments to print
#  RETURN VALUE
#    Prints arguments to $channel in C-style fprintf and flashes the $channel.
#  EXAMPLE
#    ::tclu::printfFlush $channelID "%s%d\n" "Value of i:" $i
#********
#------------------------------------------------------------------------

proc ::tclu::fprintfFlush {channel format args} {
    puts $channel [eval format [list $format] $args]
    flush $channel
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::public
#  NAME
#    ::tclu::public -- miminc a "public" keyword of [Incr Tcl]
#  USAGE
#    ::tclu::public proc1 ?proc2? ...
#
#  DESCRIPTION
#    This is dummy proc - used to enhance the syntax and
#    readability of the code. This should be used to define public procs 
#    in a namescape
#  ARGUMENTS
#    args  - a procs to be registered as public
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::public proc1 proc2 proc3
#********
#------------------------------------------------------------------------

proc ::tclu::public {args} {
    variable public 

    set public(args) $args
    uplevel 1 {	
	set namespace [namespace current]
	foreach proc $::tclu::public(args) { 
	    lappend ::tclu::public(procs,$namespace) $proc 
	}
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::getPublic
#  NAME
#    ::tclu::getPublic -- returns all public procs inside a current namespace
#                        that matches a pattern
#  USAGE
#    ::tclu::getPublic ?pattern?
#
#  DESCRIPTION
#    Returns all public procs names inside a current namespace that matches a 
#    pattern. If pattern is omitted, then all public procs manes are returned.
#  ARGUMENTS
#    pattern  - a pattern to match against the public proc names
#  RETURN VALUE
#    All public procs names inside a current namespace that matches.
#  EXAMPLE
#    ::tclu::getPublic *get*
#  SEE ALSO
#    ::tclu::public
#********
#------------------------------------------------------------------------
proc ::tclu::getPublic {{pattern *}} {
    variable public 
    
    set public(pattern) $pattern
    uplevel 1 {
	set procs {}
	set namespace [namespace current]
	foreach proc $::tclu::public(procs,$namespace) {
	    if { [string match $::tclu::public(pattern) $proc] } {
		lappend  procs $proc
	    }
	}
	return $procs
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::private
#  NAME
#    ::tclu::private -- miminc a "private" keyowrd of C++
#  USAGE
#    ::tclu::private proc1 ?proc2? ...
#
#  DESCRIPTION
#    This is dummy proc - used to enhance the syntax and
#    readability of the code. This should be used to define private procs 
#    in a namescape
#  ARGUMENTS
#    args  - a procs to be registered as private
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::private proc1 proc2 proc3
#********
#------------------------------------------------------------------------

proc ::tclu::private {args} {
    variable private 

    set private(args) $args
    uplevel 1 {	
	set namespace [namespace current]
	foreach proc $::tclu::private(args) { 
	    lappend ::tclu::private(procs,$namespace) $proc 
	}
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::getPrivate
#  NAME
#    ::tclu::getPrivate -- returns all private procs inside a current namespace
#                        that matches a pattern
#  USAGE
#    ::tclu::getPrivate ?pattern?
#
#  DESCRIPTION
#    Returns all private procs names inside a current namespace that matches a 
#    pattern. If pattern is omitted, then all private procs manes are returned.
#  ARGUMENTS
#    pattern  - a pattern to match against the private proc names
#  RETURN VALUE
#    All private procs names inside a current namespace that matches.
#  EXAMPLE
#    ::tclu::getPrivate *get*
#  SEE ALSO
#    ::tclu::private
#********
#------------------------------------------------------------------------

proc ::tclu::getPrivate {{pattern *}} {
    variable private 
    
    set private(pattern) $pattern
    uplevel 1 {
	set procs {}
	set namespace [namespace current]
	foreach proc $::tclu::private(procs,$namespace) {
	    if { [string match $::tclu::private(pattern) $proc] } {
		lappend  procs $proc
	    }
	}
	return $procs
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lsearchCase
#  NAME
#    ::tclu::lsearchCase -- case insensitive variant of Tcl lsearch command
#  USAGE
#    ::tclu::lsearchCase list pattern
#
#  DESCRIPTION
#    Case insensitive variant of Tcl lsearch command, with the exception
#    that lsearch options are omitted here.
#  ARGUMENTS
#    list    - a list to search
#    pattern - a pattern to search against
#  RETURN VALUE
#    The return value of Tcl's lsearch command.
#  EXAMPLE
#    ::tclu::lsearchCase $list *get*
#********
#------------------------------------------------------------------------

proc ::tclu::lsearchCase {list pattern} {
    return [lsearch [string toupper $list] [string toupper $pattern]]
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::range
#  NAME
#    ::tclu::range -- variation of lrange Tcl command, the difference is that it
#                     concatenates the elements
#  USAGE
#    ::tclu::range list start end
#
#  DESCRIPTION
#    It is like lrange Tcl command, the difference is that it concats 
#    (i.e. format "%s ") insted of list the elements together.
#  ARGUMENTS
#    list    - a list to search
#    start   - first element in the range
#    end     - last element in the range (can also have a value "end")
#  RETURN VALUE
#    Return the elements in the range [start,end] concatenated together.
#  EXAMPLE
#    ::tclu::range $list 1 3
#********
#------------------------------------------------------------------------

proc ::tclu::range {list start end} {
    set out {}
    if { $end == "end" } {
	set end [expr [llength $list] - 1]
    }
    for {set i $start} {$i <= $end} {incr i} {
	append out [format "%s " [lindex $list $i]]
    }
    return $out
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::newset
#  NAME
#    ::tclu::newset -- sets the variable ONLY if it doesn't already exists 
#  USAGE
#    ::tclu::newset varName value
#  DESCRIPTION
#    Sets the variable only if it doesn't already exists !!!
#  ARGUMENTS
#    varName  - name of a variable to set
#    value    - value to be assigned to varName 
#  RETURN VALUE
#    Returns $value.
#  EXAMPLE
#    ::tclu::newset today {Dec 14}
#********
#------------------------------------------------------------------------

proc ::tclu::newset { varName {value {}} } {
    upvar 1 $varName var
    if { ![info exists var] } {
	set var $value
    }
    return $value
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::conset
#  NAME
#    ::tclu::conset -- conditional setting of the variable
#  USAGE
#    ::tclu::conset var condition truevalue falseValue
#  DESCRIPTION
#     Conditional set of the variable. Similar to C construct:
#     var = condition ? trueValue : falseValue
#  ARGUMENTS
#     varName     - name of the variable
#     condition   - condition to match upon
#     trueValue   - if condition is TRUE  set var to this value
#     falseValue  - if condition is FALSE set var to this value
#  RETURN VALUE
#     Returns the value that has beeen assigned to variable.
#  EXAMPLE
#     ::tclu::conset max "$a > $b" $a $b
#  TODO
#     Maybe the following syntax should be implemented:
#     ::tclu::conset max "$a > $b" ? $a : $b --> NO, because you can
#     use: expr { textExpr ? trueValue : falseValue }
#********
#------------------------------------------------------------------------

proc ::tclu::conset {varName condition trueValue falseValue} {
    upvar 1 $varName var
    if { $condition } {
	set var $trueValue
    } else {
	set var $falseValue
    }
    return $var
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::ifincr
#  NAME
#    ::tclu::ifincr -- increases the variable by value
#  USAGE
#    ::tclu::ifincr varName ?value?
#  DESCRIPTION
#    Increase the variable by value (default for value is 1). 
#    If it doesn't already exists it sets it to zero.
#  ARGUMENTS
#    varName  - name of the variable
#    value    - increase value
#  RETURN VALUE
#    Returns the value of a variable on exit.
#  EXAMPLE
#    ::tclu::ifincr ipol -1
#********
#------------------------------------------------------------------------

proc ::tclu::ifincr {varName {value 1}} {
    upvar 1 $varName var
    if { ![info exists var] } {
	set var 0
    } else {
	incr var $value
    }
    return $var
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::ifexists
#  NAME
#    ::tclu::ifexists -- execute a script in uplevel 1 level if variable exists
#  USAGE
#    ::tclu::ifexists varName script
#  DESCRIPTION
#    Execute the script only when variable varName exists.
#  ARGUMENTS
#    varName  - name of the variable
#    script   - script to execute
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::ifexists ipol { incr ipol -2 }
#********
#------------------------------------------------------------------------

proc ::tclu::ifexists {varName script} {
    upvar 1 $varName var
    variable ifexists

    set ifexists(script) $script
    if { [info exists var] } {
	global ::tclu::ifexists
	uplevel 1 {
	    eval $::tclu::ifexists(script)
	}
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::scan
#  NAME
#    ::tclu::scan -- an enhanced scan command with the %S option
#  USAGE
#    ::tclu::scan string formatString args
#  DESCRIPTION
#    This is an enhanced scan command with an additional %S or %nS
#    option, n being an integer. The %S option means a text-string which
#    can contain the white-spaces. The %S reads up to the end of the string,
#    while %nS reads n characters.
#  ARGUMENTS
#    string        - string to scan
#    formatString  - string which holds the format conversion specifiers
#    args          - list of variables which will hold the scaned values 
#  RETURN VALUE
#    Number of scan conversions or -1 if an EOF occurs.
#  EXAMPLE
#    ::tclu::scan $string "%S %6.3f" textString floatNumber
#********
#------------------------------------------------------------------------

proc ::tclu::scan {string formatString args} {

    #------------------------------------------------------------------------
    # possible solutions for the %S
    # -----------------------------
    # 1.) %S reads the whole string or ???reads "to je tone" or 'to je tone'???
    # 2.) %[0-9]+S reads [0-9]+ characters !!!
    # hence this means one must be extremely careful with it, since
    # %s%S%f would be a nonsense, since %S would scan trough the whole string !!!
    #------------------------------------------------------------------------

    #
    # Search for the %S format-string and make "composite" scan command.
    # Example:
    #          ::tclu::scan {1.2 2.5 {To je Tone}} %15.10f%15.10f%S var1 var2 var3
    #    
    # this example should store 1.2 in var1, 2.5 in var2, and "To je
    # Tone" in var3

    set startIndex 0
    set nvar       [llength $args]

    for {set i 0} {$i < $nvar} {incr i} {
	set var  [lindex $args $i]
	upvar $var varname
	
	set _fmt [scan__NextFmtConversion $formatString]

	if { [regexp {^%[0-9]*S} $_fmt] } {
	    if { [regexp {^%[0-9]+S} $_fmt] } {
		set varname [string range $string 0 [expr [string trim $_fmt %S] - 1]]
	    } else {
		set varname $string
		return $i
	    }
	} elseif { [::scan "$string" $_fmt varname] < 0 } {
	    return -1
	}	
	set formatString [string range $formatString [string length $_fmt] end]
	set pos          [string first $varname $string]
	set string       [string range $string [expr [string length $varname] + 1 + $pos] end]
    }
    return [incr i -1]
}
proc ::tclu::scan__NextFmtConversion {string} {
    set ind [string first % $string 0]
    set ind [string first % $string [incr ind]]
    if { $ind == -1 } {
	return [string range $string 0 end]
    } else {
	return [string range $string 0 [expr $ind - 1]]
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::isElementInList
#  NAME
#    ::tclu::isElementInList -- checks if the element is in list
#  USAGE
#    ::tclu::isElementInList listVar element
#  DESCRIPTION
#    This proc checks if element is present in the list. 
#  ARGUMENTS
#    listVar  - name of the variable holding the list
#    element  - value of a list element to query
#  RETURN VALUE
#    returns 1 if element is present and 0 otherwise.
#  EXAMPLE
#    ::tclu::isElementInList thisList $element
#********
#------------------------------------------------------------------------

proc ::tclu::isElementInList {listVar element} {
    upvar $listVar listValue
    set ind [lsearch -exact $listValue $element]
    if { $ind < 0 } {
	return 0
    } else {
	return 1
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::removeElementFromList
#  NAME
#    ::tclu::removeElementFromList -- removes element from list
#  USAGE
#    ::tclu::removeElementFromList listVar element
#  DESCRIPTION
#    This proc removes the element from list if it is present.
#  ARGUMENTS
#    listVar  - name of the variable holding the list
#    element  - value of a list element to query
#  RETURN VALUE
#    returns 1 if element was removed and 0 otherwise.
#  EXAMPLE
#    ::tclu::removeElementFromList thisList $element
#********
#------------------------------------------------------------------------

proc ::tclu::removeElementFromList {listVar element} {
    upvar $listVar listValue
    set ind [lsearch -exact $listValue $element]
    if { $ind > -1 } {
	set l1 [lrange $listValue 0 [expr $ind - 1]]
	set l2 [lrange $listValue [expr $ind + 1] end]
	uplevel "set $listVar \"[concat $l1 $l2]\""
	return 1
    } else {
	return 0
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::addElementToList
#  NAME
#    ::tclu::addElementToList -- adds an element to list if it's not in the list
#  USAGE
#    ::tclu::addElementToList listVar element
#  DESCRIPTION
#    This proc adds an element to the list only when it is not already 
#    in the list.
#  ARGUMENTS
#    listVar  - name of the variable holding the list
#    element  - value of a list element to query
#  RETURN VALUE
#    returns 1 if element was added and 0 otherwise.
#  EXAMPLE
#    ::tclu::addElementToList thisList $element
#********
#------------------------------------------------------------------------

proc ::tclu::addElementToList {listVar element} {
    upvar $listVar listValue
    if { [lsearch -exact $listValue $element] < 0 } {
	uplevel "lappend $listVar [list $element]"
	return 1
    } else {
	return 0
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lindexEnd
#  NAME
#    ::tclu::lindexEnd -- like lindex, but ...
#  USAGE
#    ::tclu::lindexEnd list index
#  DESCRIPTION
# This proc is like Tcl's lindex, but if index is larger then "end"
# then it returns last element in list.
#  ARGUMENTS
#    list  - a list
#    index - index of the list element
#  RETURN VALUE
#    returns the index's element from list, or the last one if index > end.
#  EXAMPLE
#    ::tclu::lindexEnd thisList $element
#********
#------------------------------------------------------------------------

proc ::tclu::lindexEnd {list index} {
    set length [llength $list]
    set index  [expr {$index > $length ? "end" : $index}]
    return [lindex $list $index]
}

# ------------------------------------------------------------------------
#   INCOMING
# ------------------------------------------------------------------------

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::format
#  NAME
#    ::tclu::format -- an enhanced format command with the %S option 
#  USAGE
#    ::tclu::format formatString args
#  DESCRIPTION
#    This is an enhanced format command with the additional %S
# option. The %S option means a text-string which can contain the
# white-spaces. For format %s %S would be the same, but not for the
# scan. Therefore %S in format is used for compatibility of %S in
# ::tclu::scan
#  ARGUMENTS
#    formatString  - string which holds the format conversion specifiers
#    args          - arguments which will be formated
#  RETURN VALUE
#    Formated string.
#  EXAMPLE
#    ::tclu::format "%S %6.3f" {Value of the variable is:} 10.3
#********
#------------------------------------------------------------------------

proc ::tclu::format {formatString args} {
    regsub -all {(%)([-]?)([0-9]*)(S)} $formatString {\1\2\3s} formatString
    puts $formatString
    return [eval format [list $formatString] $args]
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::tempFile
#  NAME
#    ::tclu::tempFile -- utility for managing temporary files
#  USAGE
#    ::tclu::tempFile name|open|delete ?args?
#  DESCRIPTION
#    This procedure namages the temporary files. It has three different 
#  modes of operations:
#
#    1. tempFile name ?prefix? -- returns the full path-name of the temporary
#                                 file. If prefix is specified, the filename
#                                 will be of type prefix.* 
#
#    2. tempFile open ?prefix? -- does as tempFile name, but additionally opens
#                                 the file for writing. Returns the channel ID.
#    
#    3. tempFile delete all|filename ?filename1?
#                              -- deletes the temporary files. The "tempFile 
#                                 delete all" deletes all files, whereas in 
#                                 other mode only specified files are deleted.
#
#  ARGUMENTS
#    what -- the mode of operation (i.e. name, open, or delete)
#    args -- arguments according to above description
#  RETURN VALUE
#    1. In mode "name" the routine returns the full temporary-file pathname.
#    2. In mode "open" the routine returns the channel ID of opened temporary-file
#    3. In mode "delete" the routine returns an empty string.
#  EXAMPLE
#    set tempfile [::tclu::TempFile name]
#********
#------------------------------------------------------------------------

proc ::tclu::tempFile {what args} {
    variable tempFile

    switch -- $what {
	name {
	    #
	    # returns the name of a tempfile
	    #
	    # usage: tempFile name ?prefix?
	    return [_tempFile_name $args]
	}

	open {
	    #
	    # set the name and opens the tempFile for writing;
	    # returns the channel
	    #
	    # usage: tempFile open ?prefix?
	    set tmpfile [::tclu::_tempFile_name $args]
	    set channel [open $tmpfile w]
	    return $channel
	}

	delete {
	    #
	    # delete a temporary directory
	    #
	    # usage: tempFile delete all|$filenames
	    _tempFile_delete $tempFile(files) $args
	}

	default {
	    error "expected \"name\", \"open\", or \"delete\", but got \"$what\""
	}
    }
    return ""
}
proc ::tclu::_tempFile_name {args} {
    variable tempFile
    global   env
    #
    # returns the name of a tempfile
    #
    # usage: tempFile name ?prefix?
    if { $args == "" } {
	set prefix temp
    } else {
	set prefix [lindex $args 0]
    }
    
    switch $::tcl_platform(platform) {
	unix {
	    set tmpdir /tmp 
	    if { ![file writable $tmpdir] } {
		if { [file writable [pwd]] } {
		    set tmpdir [pwd]
		} else {
		    set tmpdir $env(HOME)
		}
	    }
	    set file ${prefix}.[pid].$tempFile(counter)
	} macintosh {
	    set tmpdir $env(TRASH_FOLDER) ; # a better place?
	    set file ${prefix}[pid]n$tempFile(counter)
	} default {
	    set tmpdir [pwd]	    
	    if { ! [catch {set _tmpdir $env(TMP)}] } {
		if { [file writable $_tmpdir] } { 
		    set tmpdir $env(TMP) 
		}
	    }
	    if { ! [catch {set _tmpdir $env(TEMP)}] } {
		if { [file writable $_tmpdir] } { 
		    set tmpdir $env(TEMP) 
		}
	    }
	    set file ${prefix}[pid]n$tempFile(counter)
	}
    }
    set tmpfile [file join $tmpdir $file]
    if { ![file writable $tmpdir] } {
	error "Failed to get a writable temporary directory"
    }
    
    incr tempFile(counter)
    lappend tempFile(files) $tmpfile
    return $tmpfile
}
proc ::tclu::_tempFile_delete {allfiledirs args} {
    
    #
    # delete the tempFile|tempDir
    #
    # usage: _tempFile_delete allfiledirs  all|$filenames|$dirnames
    if { $args == "" || $args == "all" } {
	set files $allfiledirs
    } else {
	set files [lindex $args 0]
    }
    catch {eval {file delete -force} $files}
    #if { $files == "all" } {
    #	# delete all tempfiles
    #	catch {eval {file delete -force} $allfiledirs}
    #} else {
    #	# delete a particular file
    #	catch {eval {file delete -force} $files}
    #}
    return ""
}

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::tempDir
#  NAME
#    ::tclu::tempDir -- utility for managing temporary directories
#  USAGE
#    ::tclu::tempDir create|delete ?args?
#  DESCRIPTION
#    This procedure namages the temporary directories. It has two different 
#  modes of operations:
#
#    1. tempDir create ?prefix? -- returns the full path-name of the created
#                                  temporary directory. If prefix is specified, 
#                                  the dirname will be of type prefix.* 
#
#    3. tempDir delete all|dirname ?dirname1?
#                              -- deletes the temporary directories. The "tempDir 
#                                 delete all" deletes all directories, whereas in 
#                                 other mode only specified directories are deleted.
#
#  ARGUMENTS
#    what -- the mode of operation (i.e. create or delete)
#    args -- arguments according to above description
#  RETURN VALUE
#    1. In mode "create" the routine returns the full temporary-dir pathname.
#    3. In mode "delete" the routine returns an empty string.
#  EXAMPLE
#    set tempfile [::tclu::TempDir name]
#********
#------------------------------------------------------------------------

proc ::tclu::tempDir {what args} {
    variable tempFile

    switch -- $what {
	create {
	    #
	    # create a temporary directory
	    #
	    # usage: tempDir create ?prefix?
	    set tempdir [_tempFile_name $args]
	    file mkdir $tempdir
	    lappend tempFile(dirs) $tempdir
	    return $tempdir
	}

	delete {
	    #
	    # delete a temporary directory
	    #
	    # usage: tempDir delete all|dirname
	    _tempFile_delete $tempFile(dirs) $args 
	}

	default {
	    error "expected \"name\", \"open\", or \"delete\", but got \"$what\""
	}
    }
    return ""
}	    


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::stringMatch
#  NAME
#    ::tclu::stringMatch -- "string match" according to nocase specification
#  USAGE
#    stringMatch pattern string ?nocase?
#  DESCRIPTION
#    This procedure is a variation of the theme on Tcl's "string match
# ?-nosace? pattern string" command. The only difference between the two is
# the usage. Namely:
# 
# Tcl command:    string match ?-nocase? pattern string
# 
# ::tclu:: proc:  stringMatch pattern string 0|1
#
#  ARGUMENTS
#    pattern --
#    string  --
#    nocase  -- flag for case sensitive/insensitive match (must be 0|1)
#  RETURN VALUE
#    Returns 1 if pattern and string matches, and 0 otherwise.
#  EXAMPLE
#    set match [::tclu::stringMatch what* WhatEver 0]
#********
#------------------------------------------------------------------------

proc ::tclu::stringMatch {pattern string {nocase 0}} {
    if { $nocase } {
	return [string match -nocase $pattern $string]
    } else {
	return [string match $pattern $string]
    }
}
 

#------------------------------------------------------------------------
#****f* ::tclu/::tclu::nonblocking
#  NAME
#    ::tclu::nonblocking -- executes a program in nonblocking mode
#  USAGE
#    ::tclu::nonblocking open|exec|unset|kill ?args?
#  DESCRIPTION
#    This procedure namages the executions of external programs in
# nonblocking mode.  modes of operations:
# 
#    1. nonblocking open          -- returns an ID for subsequent call to 
#                                    "nonblocking exec"
# 
#    2. nonblocking exec id args  -- executes external program ($args) in 
#                                    nonblocking mode. The "id" is the 
#                                    return-value of the previous 
#                                    "nonblocking open" call.
#    
#    3. nonblocking unset id      -- will unset all nonblocking(*,$id) array 
#                                    elements. Call this when a nonblocking id'th
#                                    information is not needed anymore.
#
#    4. nonblocking kill id       -- kills the nonblocking id process
#
#  RETURN VALUE
#    1. In mode "open" the routine returns the ID.
#  EXAMPLE
#    set id [::tclu::nonblocking open]
#    ::tclu::nonblocking exec $myProgram $myProgramOption
#    # here I do something with the ::tclu::nonblocking(*,$id) information
#    ....
#    # now the ::tclu::nonblocking(*,$id) information is not needed anymore
#    ::tclu::nonblocking unset $id
#********
#------------------------------------------------------------------------

proc ::tclu::nonblocking {mode args} {
    variable nonblocking

    switch -exact -- $mode {
	open {
	    return [incr nonblocking(counter)]
	}

	exec {
	    set count [lindex $args 0]
	    set args  [lrange $args 1 end]
	    if { [info exists nonblocking(done,$count)] } {
		::tclu::ERROR "statement \"::tclu::nonblocking exec\" executed before \"::tclu::nonblocking open\""
		return
	    }
	    set nonblocking(done,$count) 0
	    set nonblocking(fID,$count)  [open [concat | $args] r]
	    fconfigure $nonblocking(fID,$count) -blocking 0
	    fileevent  $nonblocking(fID,$count) readable [list ::tclu::_nonblockingEvent $count]
	    
	    tkwait variable ::tclu::nonblocking(done,$count)    
	    #::tclu::DEBUG nonblocking: releasing tkwait variable nonblocking(done,$count) ...
	    return $count
	}

	unset {
	    foreach id $args {
		if { ![info exists nonblocking(done,$id)] } {
		    ::tclu::ERROR "nonblocking exec id $id does not exists"
		}		
		array unset nonblocking *,$id
	    }
	}

	kill {
	    foreach id $args {
		if { ![info exists nonblocking(done,$id)] } {
		    ::tclu::ERROR "nonblocking exec id $id does not exists"
		}
		catch {close $nonblocking(fID,$count)}
	    }
	}
    }
}
proc ::tclu::_nonblockingEvent {count} {
    variable nonblocking
    #::tclu::DEBUG nonblocking: Event ...
    if { ! [eof $nonblocking(fID,$count)] } {
	append nonblocking(output,$count) [gets $nonblocking(fID,$count)]\n
    } else {
	if { [catch {close $nonblocking(fID,$count)}] } {
	    ::tclu::errorDialog "error while closing nonblocking execution \#.$count"
	}
	set nonblocking(done,$count) 1
    }
}


# 
# proc ::tclu::format {formatString args} {
#     
#     #
#     # Search for the %S string and make "composite" format command.
#     # Example:
#     #          ::tclu::format %15.10f%15.10f%S%s $arg1 $arg2 $arg3 $arg4
#     #
#     # will be tranformed to:
#     #
#     #          set fmtString    [::format %15.10f%15.10f $arg1 $arg2]
#     #          append fmtString $arg3
#     #          append fmtString [::format %s $arg4]
#     #
# 
#     set _old_ind   0
#     set startIndex 0
#     set strLength  [string length $formatString]
#     set output     {}
# 
#     while {1} {
# 	set ind [string first %S $formatString $startIndex]	
# 	if { $ind == -1 } {
# 	    set _fmt    [string range $formatString $startIndex end]
# 	    set _fmtStr [lrange $args $_old_ind end]
# 	    append output [eval {::format $_fmt} $_fmtStr]
# 	    break
# 	}
# 	
# 	set _fmt        [string range $formatString $startIndex [expr $ind - 1]]
# 	set _nofmt_ind  [format__whichArg $formatString $ind] 
# 	set _fmtArgs    [lrange $args $_old_ind [expr $_nofmt_ind - 1]]
# 	set _fmtStr     [eval {::format $_fmt} $_fmtArgs]
# 	set _nofmtStr   [lindex $args $_nofmt_ind]
# 	
# 	append output  ${_fmtStr}${_nofmtStr}
# 
# 	set startIndex [expr $ind + 2]
# 	set _old_ind   [expr $_nofmt_ind + 1]
#     }
# 
#     return $output
# }
# proc ::tclu::format__whichArg {formatString endIndex} {
#     set fmtStr [string range $formatString 0 $endIndex]
#     set which -1
#     set startIndex 0
#     while {1} {
# 	set ind [string first % $fmtStr $startIndex]
# 	if { $ind == -1 } {
# 	    break
# 	}
# 	incr ind
# 	if { [string range $fmtStr $ind $ind] != "%" } {
# 	    incr which
# 	}
# 	set startIndex $ind
#     }
#     return $which
# }
