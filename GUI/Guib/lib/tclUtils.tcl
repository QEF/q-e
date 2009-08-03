#
# $RCSfile: tclUtils.tcl,v $ --
#
#      This file contains the Tone Kokalj's Tcl utilities functions.
#
# Copyright (c) 2003-2008  Anton Kokalj   Email: tone.kokalj@ijs.si
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
# $Id: tclUtils.tcl,v 1.16 2009-08-03 14:10:18 kokalj Exp $ 
#

#------------------------------------------------------------------------
#****h* TclTkLib/::tclu
#  NAME
#    TCLU == Tone Kokalj's Tcl Utilities 
#                          
#  COPYRIGHT
#    2001--2004 (c) by Tone Kokalj
#  AUTHOR
#    Tone Kokalj
#  CREATION DATE
#    Starting on Fri Dec 14 08:32:59 CET 2001
#  NOTES
#    Here are the basic Tcl-only utility functions.
#******
#------------------------------------------------------------------------

package require fileutil
package provide tclu 0.9

namespace eval ::tclu {
    variable DEBUG      0
    variable DEBUG_FILE 0
    variable debug
    variable error
    variable public  
    variable private 
    variable _line "------------------------------------------------------------------------"
    variable tempFile
    variable nonblocking
    
    set tempFile(counter) 0
    set tempFile(files)   {}
    set tempFile(dirs)    {}
    set nonblocking(counter) -1

    namespace export dummy
    namespace export absolutePath
    namespace export usage
    namespace export writeFile
    namespace export readFile
    namespace export lineread
    namespace export DEBUG
    namespace export ERROR
    namespace export errorDialog
    namespace export warningDialog
    namespace export labelMsg
    namespace export putsFlush
    namespace export printf
    namespace export printfFlush
    namespace export fprintf
    namespace export fprintfFlush
    namespace export public
    namespace export getPublic
    namespace export private
    namespace export getPrivate
    namespace export range
    namespace export newset
    namespace export ifexists
    namespace export lpresent
    namespace export lremove
    namespace export lpop
    namespace export ladd
    namespace export lget
    namespace export lequal
    namespace export tempFile
    namespace export tempDir
    namespace export stringMatch
    namespace export nonblocking
    namespace export expandArgs
    namespace export getArgs
    namespace export extractArgs
    namespace export modePatternArgs
    namespace export mustBeEvenNumOfArgs
    namespace export mustBeOddNumOfArgs

    # Don't export the following cmd's, because their names clash with
    # Tcl-commands:
    #   scan
    #   format
    #   incr
}


#****f* ::tclu/::tclu::dummy
# SYNOPSIS
proc ::tclu::dummy {} {
    # PURPOSE
    #   Do nothing.
    # SOURCE
    return ""
}
#******


#****f* ::tclu/::tclu::absolutePath
# SYNOPSIS
proc ::tclu::absolutePath {path} {
    # PURPOSE
    #   Return absolute path. Please use "file normalize instead".
    # ARGUMENTS
    # * name -- name of the file
    # SOURCE
    return [file normalize $path]
}
#******


#****f* ::tclu/::tclu::usage
#  SYNOPSIS
proc ::tclu::usage {args} {
    # USAGE
    #   ::tclu::usage ?channelID? text
    # PURPOSE
    #   Prints a USAGE message. 
    # ARGUMENTS
    # * channelID -- (optional) a channel identifier as returned from call to open or socket
    # * text -- a usage-text to be printed
    # EXAMPLE
    #   ::tclu::usage stderr "::tclu::usage ?channelID? text"
    # SOURCE

    variable _line

    if { [llength $args] > 2 } {
	::tclu::usage "::tclu::usage ?channelID? text"
    }
    set channel stderr
    set text [lindex $args end]
    if { [llength $args] == 2 } {
	set channel [lindex $args 0]
    }

    puts $channel $_line
    puts $channel [labelMsg Usage $text]
    putsFlush $channel $_line\n
}
#********


#****f* ::tclu/::tclu::writeFile
# SYNOPSIS
proc ::tclu::writeFile {filename content {flag w}} {
    #  DESCRIPTION
    #    Write a content to a file.
    #  ARGUMENTS
    #  * filename  -- name of the file to be written
    #  * content   -- content to be writtento a file
    #  * flag      -- (optional) access flag for "open" Tcl command, i.e., w, w+, a, a+
    #  EXAMPLE
    #    ::tclu::writeFile test.txt {This is a test.txt file} w
    #  SOURCE 
    set   fID  [open $filename $flag]
    puts  $fID $content
    flush $fID
    close $fID
}
#********


#****f* ::tclu/::tclu::readFile
# SYNOPSIS
proc ::tclu::readFile {args} {
    #  USAGE
    #    ::tclu::readFile ?-nonewline? filename
    #
    #  DESCRIPTION
    #    Reads a file and returns its content. If the -nonewline
    #    switch is specified then the last character of the file is
    #    discarded if it is a newline.
    #
    #  SOURCE

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
#********


#****f* ::tclu/::tclu::lineread
# SYNOPSIS
proc ::tclu::lineread {var file script} {
    # PURPOSE
    #   Read entire file line-by-line and at each line execute a
    #   script at one level up.
    # ARGUMENTS
    # * var    -- name of variable where the content of line will be stored
    # * file   -- name of file to read
    # * script -- script to execute when line is read 
    #
    # CREDITS
    #   Based on fileutils::foreachLine from tcllib (almost verbatim).
    # SOURCE
    upvar $var line

    set fid    [open $file r]
    set code   0
    set result {}

    while { ! [eof $fid] } {
        gets $fid line
        set code [catch {uplevel 1 $script} result]
        if {($code != 0) && ($code != 4)} { 
            break 
        }
    }
    close $fid

    if { ($code == 0) || ($code == 3) || ($code == 4) } {
        return $result
    }
    if { $code == 1 } {
        global errorCode errorInfo
        return \
            -code      $code      \
            -errorcode $errorCode \
            -errorinfo $errorInfo \
            $result
    }
    return -code $code $result
}
#******


#****f* ::tclu/::tclu::DEBUG
# SYNOPSIS
proc ::tclu::DEBUG {args} {
    #  USAGE
    #    ::tclu::DEBUG ?-stderr? message
    #
    #  DESCRIPTION
    #    Prints a DEBUG message either to stdout or to stderr (when
    #    -stderr option is specified).
    # 
    #    If ::tclu::DEBUG_FILE variable is set it also prints debugging
    #    message to DEBUG file.
    #
    #  EXAMPLE
    #    ::tclu::DEBUG -stderr "Variable A has a value: $A"
    #  SOURCE

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
#********


#****f* ::tclu/::tclu::ERROR
# SYNOPSIS
proc ::tclu::ERROR {errMsg {info {}} {code {}}} {    
    #  DESCRIPTION
    #    Like a Tcl-error cmd, but add the ERROR: prefix to the error message
    #  ARGUMENTS
    #  * errMsg -- error message
    #  * info -- used to initialize the global variable errorInfo,
    #    which is used to accumulate a stack trace
    #  * error -- its value is stored in Tcl's errorCode variable
    #  RETURN VALUE
    #    Generates and error message.
    #  EXAMPLE
    #    ::tclu::ERROR {an error hac occurred because variable xyz is not defined}
    #  SOURCE

    variable error
    variable _line
    set error(errMsg) "[labelMsg ERROR $errMsg]\n"
    set error(info) $info
    set error(code) $code
    uplevel 1 {
	error $::tclu::error(errMsg) $::tclu::error(info) $::tclu::error(code)
    }
}
#********


#****f* ::tclu/::tclu::abort
# SYNOPSIS
proc ::tclu::abort {msg {exit_status 1}} {
    # PURPOSE   
    #   Print a message and abort (i.e. exit from application).
    # ARGUMENTS   
    # * msg -- message to print to stderr before exiting.
    # * exit_status -- (optional, default = 1) 
    # SOURCE
    puts stderr "\n[labelMsg ABORT $msg]"

    if { [info level] > 1 } {
        puts stderr "\nABORT statement executed in: [info level -1]"
    } else {
        puts stderr "\nABORT statement executed in top-level"
    }
    exit $exit_status
}


#******

#****f* ::tclu/::tclu::errorDialog
# SYNOPSIS
proc ::tclu::errorDialog {errMsg} {
    # DESCRIPTION
    #   Prints an error message either to stderr or into tk_messageBox
    #   when Tk package is present.
    # ARGUMENTS
    # * errMsg -- error message
    # EXAMPLE
    #   ::tclu::errorDialog "an error hac occurred because variable xyz is not define"

    _error_or_warning_dialog error $errMsg
}
#********


#****f* ::tclu/::tclu::warningDialog
# SYNOPSIS
proc ::tclu::warningDialog {warnMsg} {
# DESCRIPTION
#   Prints a warning message either to stderr or into tk_messageBox
#   when Tk package is present.
# ARGUMENTS
# * warnMsg -- warning message
# EXAMPLE
#   ::tclu::warningDialog "file \"$file\" does not exist"
    
    _error_or_warning_dialog warning $warnMsg
}
#********


proc ::tclu::_error_or_warning_dialog {type msg} {
    variable _line
    
    set types {error warning}

    if { [lsearch -exact $types $type] >= 0 } {

	set TYPE [string toupper $type]	
	
	if { [catch {package present Tk}] } {
	    # make a nicely formated message
	    set message    "\n$_line\n"
	    append message "[labelMsg $TYPE $msg]\n"
	    append message "$_line\n"
	    
	    puts stderr $message
	    flush stderr
	} else {
	    tk_messageBox -title $TYPE -message $msg -type ok -icon $type
	}
    } else {
	error "wrong type \"$type\" in ::tclu::_error_or_warning_dialog, should be [join $types ,]"
    }
}


#****f* ::tclu/::tclu::labelMsg
# SYNOPSIS
proc ::tclu::labelMsg {label msg} {
    # PURPOSE
    #   Create a nicely formated (multi-line) message with starting
    #   "label". In particular, the 2nd and subsequent lines will be
    #   properly indented. See the example below.
    #
    # EXAMPLE
    # This call:
    #   ::tclu::labelMsg ERROR { the soure of the error is
    #        a bit curious. 
    #   }
    #
    # will produce the following message:
    #
    #    ERROR: the source of the error is
    #           a bit curious.
    #
    # ARGUMENTS    
    # * label -- label-text to apper as a prefix of the message
    # * msg   -- message to label and format    
    # SOURCE

    set il 1
    set len [string length $label]
    set message {}
    foreach line [split [string trim $msg] \n] {
	if { $il == 1 } {
	    append message [::format "%${len}s: %s\n" $label $line]
	    incr len 
	    incr il
	} else {
	    append message [::format "%${len}s %s\n" {} $line]
	}
    }
    return [string trim $message]
}
#******


#****f* ::tclu/::tclu::putsFlush
#  SYNPOSIS			
proc ::tclu::putsFlush {args} {
    #  USAGE
    #    ::tclu::putsFlush ?-nonewline? ?channelId? string
    #
    #  DESCRIPTION
    #    Identical to Tcl's puts, but invoke the flush immediately after.
    #    See puts man-page of Tcl.
    #
    #  SOURCE

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
#********


#****f* ::tclu/::tclu::printf
#  SYNOPSIS
proc ::tclu::printf {format args} {
    #  USAGE
    #    ::tclu::printf format ?arg? ?arg? ...
    #
    #  DESCRIPTION
    #    A C-style printf: it prints to stdout like the C-style printf.
    #  ARGUMENTS
    #  *  format  -- format description
    #  *  args    -- arguments to print
    #  EXAMPLE
    #    ::tclu::printf "%s%d\n" "Value of i:" $i
    #  SOURCE

    puts stdout [eval format [list $format] $args]
}
#********


#****f* ::tclu/::tclu::printfFlush
#  SYNOPSIS
proc ::tclu::printfFlush {format args} {
    #  USAGE
    #    ::tclu::printfFlush format ?arg? ?arg? ...
    #
    #  DESCRIPTION
    #    Prints to stdout in like C-style printf and flushes stdout.
    #
    #  ARGUMENTS
    #  *  format -- format description
    #  *  args   -- arguments to print
    #  EXAMPLE
    #    ::tclu::printfFlush "%s%d\n" "Value of i:" $i
    #  SOURCE
    
    puts stdout [eval format [list $format] $args]
    flush stdout
}
#********


#****f* ::tclu/::tclu::fprintf
# SYNOPSIS
proc ::tclu::fprintf {channel format args} {
    #  USAGE
    #    ::tclu::fprintf channel format ?arg? ?arg? ...
    #  DESCRIPTION
    #    A C-style fprintf: prints to $channel like C-style fprintf would.
    #  ARGUMENTS
    #  * channel -- a channel identifier such as returned from call to 
    #               open or socket
    #  * format  -- format description
    #  * args    -- arguments to print
    #  EXAMPLE
    #    ::tclu::fprintf $channelID "%s%d\n" "Value of i:" $i
    #  SOURCE

    puts $channel [eval format [list $format] $args]
}
#********


#****f* ::tclu/::tclu::fprintfFlush
# SYNOPSIS
proc ::tclu::fprintfFlush {channel format args} {
    #  USAGE
    #    ::tclu::fprintf channel format ?arg? ?arg? ...
    #  DESCRIPTION
    #    Prints to $channel like the C-style fprintf would and flashes the $channel.
    #  ARGUMENTS
    #  * channel -- a channel identifier such as returned from call to 
    #               open or socket
    #  * format  -- format description
    #  *  args   -- arguments to print
    #  EXAMPLE
    #    ::tclu::printfFlush $channelID "%s%d\n" "Value of i:" $i
    #  SOURCE
    
    puts $channel [eval format [list $format] $args]
    flush $channel
}
#********


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::public
#  NAME
#    ::tclu::public -- miminc a "public" keyword of [Incr Tcl]
#  USAGE
#    ::tclu::public proc1 ?proc2? ...
#
#  DESCRIPTION
#    This is dummy proc, used to enhance the syntax and
#    readability of the code. This should be used to define public procs 
#    in a namescape.
#  ARGUMENTS
#  * args -- a procs to be registered as public
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
#                         that matches a pattern
#  USAGE
#    ::tclu::getPublic ?pattern?
#
#  DESCRIPTION
#    Returns all public procs names inside a current namespace that matches a 
#    pattern. If pattern is omitted, then all public procs manes are returned.
#  ARGUMENTS
#  *  pattern -- a pattern to match against the public proc names
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
#  * args -- a procs to be registered as private
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
#  * pattern  -- a pattern to match against the private proc names
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
#  * list   -- a list to search
#  * start  -- first element in the range
#  * end    -- last element in the range (can also have a value "end")
#  RETURN VALUE
#    Return the elements in the range [start,end] concatenated together.
#  EXAMPLE
#    ::tclu::range $list 1 3
#********
#------------------------------------------------------------------------

proc ::tclu::range {list start end} {
    set out {}
    if { $end == "end" } {
	set end [expr {[llength $list] - 1}]
    }
    for {set i $start} {$i <= $end} {incr i} {
	append out [format "%s " [lindex $list $i]]
    }
    return [string trim $out]
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
#****f* ::tclu/::tclu::incr
#  NAME
#    ::tclu::incr -- increases the variable by value
#  USAGE
#    ::tclu::incr varName ?value?
#  DESCRIPTION
#    Like the Tcl incr command, but if the variable doesn't already
#    exists it sets it to zero.
#  ARGUMENTS
#  * varName  -- name of the variable
#  * value    -- increase value
#  RETURN VALUE
#    Returns the value of a variable on exit.
#  EXAMPLE
#    ::tclu::incr ipol -1
#********
#------------------------------------------------------------------------

proc ::tclu::incr {varName {value 1}} {
    upvar 1 $varName var
    if { ![info exists var] } {
	set var 0
    } else {
	::incr var $value
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
#  *  varName  -- name of the variable
#  *  script   -- script to execute
#  RETURN VALUE
#    None.
#  EXAMPLE
#    ::tclu::ifexists ipol { incr ipol -2 }
#********
#------------------------------------------------------------------------

proc ::tclu::ifexists {varName script} {
    upvar 1 $varName var

    if { [info exists var] } {
	uplevel 1 $script
    }    
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::scan
#  NAME
#    ::tclu::scan -- an enhanced Tcl scan command with the %S option
#  USAGE
#    ::tclu::scan string formatString var1 ?var2? ...
#  DESCRIPTION
#    This is an enhanced Tcl scan command with an additional %S or %nS
#    formating option, n being an integer. The %S option means a text-string which
#    can contain the white-spaces. The %S reads up to the end of the string,
#    while %nS reads n characters.
#  ARGUMENTS
#  *  string       -- string to scan
#  *  formatString -- string which holds the format conversion specifiers
#  *  args         -- a list of variables which will hold the scaned values 
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
	set varname  [lindex $args $i]
	upvar $varname var
	
	set _fmt [_nextFmtSpec $formatString]

	if { [regexp {^%[0-9]*S} $_fmt] } {

	    # %S formating option has been specified

	    if { [regexp {^%[0-9]+S} $_fmt] } {
		# scan only requested number of characters, and assing them to var
		set var [string range $string 0 [expr {[string trim $_fmt %S] - 1}]]
	    } else {
		# scan the whole string, i.e. assign $string to var
		set var $string
		return [expr {$i + 1}]
	    }
	} elseif { [::scan "$string" $_fmt var] < 0 } {
	    return -1
	}	
	set formatString [string range $formatString [string length $_fmt] end]
	set pos          [string first $var $string]
	set string       [string range $string [expr {[string length $var] + 1 + $pos}] end]
    }
    return $i
}

proc ::tclu::_nextFmtSpec {string} {
    set ind [string first % $string 0]
    set ind [string first % $string [incr ind]]
    if { $ind == -1 } {
	return [string range $string 0 end]
    } else {
	return [string range $string 0 [expr {$ind - 1}]]
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lpresent
#  NAME
#    ::tclu::lpresent -- checks if the element is in list
#  USAGE
#    ::tclu::lpresent list element
#  DESCRIPTION
#    This proc checks if element is present in the list. 
#  ARGUMENTS
#  * list     -- list to check
#  * element  -- value of a list element to query
#  RETURN VALUE
#    Returns 1 if element is present in list and 0 otherwise.
#  EXAMPLE
#    ::tclu::lpresent $list $element
#********
#------------------------------------------------------------------------

proc ::tclu::lpresent {list element} {
    set ind [lsearch -exact $list $element]
    return [expr {$ind >= 0 ? 1 : 0}]
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lremove
#  NAME
#    ::tclu::lremove -- removes element from list
#  USAGE
#    ::tclu::lremove listVar element
#  DESCRIPTION
#    This proc removes all occurrences on the element from the listVar.
#  ARGUMENTS
#  * listVar  -- name of the variable holding the list
#  * element  -- value of a list element to query
#  RETURN VALUE
#    returns 1 if element was removed and 0 otherwise.
#  EXAMPLE
#    ::tclu::lremove thisList $element
#********
#------------------------------------------------------------------------

proc ::tclu::lremove {listVar element} {
    upvar $listVar lvar

    # Fast-check: does listVar holds a single element which is equal to $element ?

    if { $lvar == $element } {
	# remove the only element
	set lvar {}
	return 1
    }
    
    set result 0
    foreach elem $lvar {
	if { ! [string match $element $elem] } {
	    lappend nl $elem
	} else {
	    set result 1
	}
    }
    if { [info exists nl] } {
	set lvar $nl
    }

    return $result
}



#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lpop
#  NAME
#    ::tclu::lpop -- remove i-th element from list
#  USAGE
#    ::tclu::lpop listVar index
#  DESCRIPTION
#    This proc removes the ith element from the listVar.
#  ARGUMENTS
#  * listVar  -- name of the variable holding the list
#  * index    -- index (0 ... end) of element to drop from list [optional, default = "end"]
#  RETURN VALUE
#    New list with removed element.
#  EXAMPLE
#    ::tclu::lpop thisList end
#********
#------------------------------------------------------------------------

proc ::tclu::lpop {listVar {index end}} {
    upvar $listVar lvar
    
    set len [llength $lvar]
    
    if { $index == "end" } {
	set index [expr {$len - 1}]
    }
    
    if { ! [string is integer $index] } {
	error "wrong index $index, must be integer number or \"end\""
    }
    
    if { $index < 0 || $index > $len - 1 } {
	error "index out of range"
    }
    
    set result ""
    set count  0
    
    foreach elem $lvar {	
	if { $count != $index } {
	    lappend result $elem
	}
	incr count
    }
    set lvar $result    
    
    return $result
}



#------------------------------------------------------------------------
#****f* ::tclu/::tclu::ladd
#  NAME
#    ::tclu::ladd -- adds an element to list only if it's not in the list
#  USAGE
#    ::tclu::ladd listVar element
#  DESCRIPTION
#    This proc is similar to Tcl lappend, put it adds an element to
#    the list only when it is not already present the list.
#  ARGUMENTS
#  * listVar -- name of the variable holding the list
#  * element -- value of a list element to query
#  RETURN VALUE
#    returns 1 if element was added and 0 otherwise.
#  EXAMPLE
#    ::tclu::ladd thisList $element
#********
#------------------------------------------------------------------------

proc ::tclu::ladd {listVar element} {
    upvar $listVar lvar
    if { [lsearch -exact $lvar $element] < 0 } {
	lappend lvar $element
	#uplevel "lappend $listVar [list $element]"
	return 1
    } else {
	return 0
    }
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lget
#  USAGE
#    ::tclu::lget list index
#  DESCRIPTION
#   This proc is like Tcl lindex, but if index is larger then "end"
#   then it returns last element in list.
#  ARGUMENTS
#  * list  -- a list
#  * index -- index of the list element
#  RETURN VALUE
#    returns the index's element from list, or the last one if index > end.
#  EXAMPLE
#    ::tclu::lget $mylist $element
#********
#------------------------------------------------------------------------

proc ::tclu::lget {list index} {
    set length [llength $list]
    set lm     [expr {$length - 1}]
    set index  [expr {$index > $lm ? "end" : $index}]
    return [lindex $list $index]
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::lequal
#  USAGE
#    ::tclu::lequal list1 list2
#  DESCRIPTION
#   This proc checks if list1 is equal to list2 and if so it returns
#   1, otherwise 0.
#  ARGUMENTS
#  * list1  -- 1st list
#  * list2  -- 2nd list
#  RETURN VALUE
#    1 or 0, depending on the equality of list 1 and list2.
#********
#------------------------------------------------------------------------

proc ::tclu::lequal {list1 list2} {
    if { [llength $list1] != [llength $list2] } {
	return 0
    }
    foreach e1 $list1 e2 $list2 {
	if { ! [string match $e1 $e2] } { return 0 }
    }
    return 1
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::format
#  NAME
#    ::tclu::format -- an enhanced Tcl-format command with the %S option added
#  USAGE
#    ::tclu::format formatString args
#  DESCRIPTION
#    This is an enhanced format command with the additional %S
#    option. The %S option means a text-string which can contain the
#    white-spaces. For format %s and %S would be the same, but not for the
#    scan. Therefore %S in format is used for compatibility of %S in
#    ::tclu::scan
#  ARGUMENTS
#  * formatString -- string which holds the format conversion specifiers
#  * args         -- arguments which will be formated
#  RETURN VALUE
#    Formated string.
#  EXAMPLE
#    ::tclu::format "%S %6.3f" {Value of the variable is:} 10.3
#********
#------------------------------------------------------------------------

proc ::tclu::format {formatString args} {
    regsub -all {(%)([-]?)([0-9]*)(S)} $formatString {\1\2\3s} formatString
    return [eval ::format [list $formatString] $args]
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::tempFile
#  NAME
#    ::tclu::tempFile -- utility for managing temporary files
#  USAGE
#    ::tclu::tempFile name|open|delete ?args?
#  DESCRIPTION
#    This procedure namages the temporary files. It has three different 
#    modes of operations:
#
#       tempFile name ?prefix? -- returns the full path-name of the temporary
#                                 file. If prefix is specified, the filename
#                                 will be of type prefix.* 
#
#       tempFile open ?prefix? -- does as tempFile name, but additionally opens
#                                 the file for writing. Returns the channel ID.
#    
#       tempFile delete all
#       or
#       tempFile delete $filename1 ?$filename2?
#                              -- deletes the temporary files. The "tempFile 
#                                 delete all" deletes all temporary files, whereas in 
#                                 other mode only specified files are deleted.
#
#  RETURN VALUE
#    1. In mode "name" the routine returns the full temporary-file pathname.
#    2. In mode "open" the routine returns the channel ID of opened temporary-file
#    3. In mode "delete" the routine returns an empty string.
#  EXAMPLE
#    set tempfile [::tclu::tempFile name]
#********
#------------------------------------------------------------------------

proc ::tclu::tempFile {what args} {
    variable tempFile

    switch -- $what {
	name {
	    #
	    # returns the name of a tempfile
	    #
	    # Usage: tempFile name ?prefix?
	    return [_tempFile_name $args]
	}

	open {
	    #
	    # set the name and opens the tempFile for writing;
	    # returns the channel
	    #
	    # Usage: tempFile open ?prefix?
	    set tmpfile [::tclu::_tempFile_name $args]
	    set channel [open $tmpfile w]
	    return $channel
	}

	delete {
	    #
	    # delete a temporary directory
	    #
	    # Usage: tempFile delete all|$filenames
	    if { $args == "" } {
		::tclu::usage "::tclu::tempFile delete all\nor\n::tclu::tempFile delete \$filename1 ?\$filename2? ..."
	    }
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
    global   env tcl_platform
    # PURPOSE
    #   Returns the name of a tempfile.    
    # USAGE
    #   _tempFile_name ?prefix?
    if { $args == "" } {
	set prefix temp
    } else {
	set prefix [lindex $args 0]
    }
    
    set tmpdir [::fileutil::tempdir]
    if { ! [file writable $tmpdir] } {
	error "Failed to get a writable temporary directory"
    }
    
    if { $tcl_platform(platform) == "unix" } {
	set file ${prefix}.[pid].$tempFile(counter)
    } else {
	set file ${prefix}[pid]n$tempFile(counter)
    }

    set tmpfile [file join $tmpdir $file]
    
    incr tempFile(counter)
    lappend tempFile(files) $tmpfile

    return $tmpfile
}

proc ::tclu::_tempFile_delete {allfiledirs what} {
    
    #
    # delete the tempFile|tempDir
    #
    # usage: _tempFile_delete allfiledirs  all|$filenames|$dirnames
    
    if { $what == "all" } {
	set files $allfiledirs
    } else {
	set files $what
    }
    catch {eval {file delete -force} $files}
}


#------------------------------------------------------------------------
#****f* ::tclu/::tclu::tempDir
#  NAME
#    ::tclu::tempDir -- utility for managing temporary directories
#  USAGE
#    ::tclu::tempDir create|delete ?args?
#  DESCRIPTION
#  This procedure namages the temporary directories. It has two different 
#  modes of operations:
#
#       tempDir create ?prefix? -- returns the full path-name of the created
#                                  temporary directory. If prefix is specified, 
#                                  the dirname will be of type prefix.* 
#
#       tempDir delete all
#       or
#       tempDir delete $dirname1 ?$dirname2? ...
#                              -- deletes the temporary directories. The "tempDir 
#                                 delete all" deletes all directories, whereas in 
#                                 other mode only specified directories are deleted.
#
#  RETURN VALUE
#  * In mode "create" the routine returns the full temporary-dir pathname.
#  * In mode "delete" the routine returns an empty string.
#  EXAMPLE
#    set tmpdir [::tclu::TempDir $name]
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
	    if { $args == "" } {
		::tclu::usage "::tclu::tempDir delete all\nor\n::tclu::tempDir delete \$dirname1 ?\$dirname2? ..."
	    }
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
#    ?-nosace? pattern string" command. The only difference between the two is
#    the usage. Namely:
# 
# Tcl command:    string match ?-nocase? pattern string
# 
# ::tclu:: proc:  stringMatch pattern string 0|1
#
#  ARGUMENTS
#  *  pattern -- as in "string match $pattern $string"
#  *  string  -- as in "string match $pattern $string"
#  *  nocase  -- flag for case sensitive/insensitive match (must be 0|1)
#  RETURN VALUE
#    Returns 1 if pattern and string matches, and 0 otherwise.
#  EXAMPLE
#    set match [::tclu::stringMatch what* WhatEver 0]
#********
#------------------------------------------------------------------------

proc ::tclu::stringMatch {pattern string {nocase 0}} {
    if { $nocase == 0 } {
	return [string match $pattern $string]
    } else {
	return [string match -nocase $pattern $string]
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
#    nonblocking mode. Modes of operations:
# 
#    1. nonblocking open          -- returns an ID for subsequent call to "nonblocking exec"
# 
#    2. nonblocking exec id args  -- executes external program ($args) in nonblocking mode. The "id" is the return-value of the previous "nonblocking open" call.
#    
#    3. nonblocking stdout id cmd -- will call the cmd and pass the new text on stdout to a proc "cmd". The proc should be of form:  cmd id text
#
#    4. nonblocking save id file  -- will save the stdout to file
#
#    5. nonblocking unset id      -- will unset all nonblocking(*,$id) array elements. Call this when a nonblocking id'th information is not needed anymore.
#
#    6. nonblocking kill id       -- kills the nonblocking id process
#
#  RETURN VALUE
#    1. In mode "open" the routine returns the ID.
#
#  EXAMPLE
#  Here is a simple example:
#    set id [::tclu::nonblocking open]
#    ::tclu::nonblocking exec $id $myProgram $myProgramOption
#
#    # here do something with the ::tclu::nonblocking(*,$id) information
#    ....
#
#    # now the ::tclu::nonblocking(*,$id) information is not needed anymore
#    ::tclu::nonblocking unset $id
#********
#------------------------------------------------------------------------

proc ::tclu::nonblocking {mode args} {
    variable nonblocking

    switch -exact -- $mode {
	open {
	    set count [incr nonblocking(counter)]
	    set nonblocking(status,$count) 1
	    return $count
	}

	exec {
	    set count [lindex $args 0]
	    if { ! [string is integer $count] } {
		::tclu::ERROR "expected integer but got $count"
		return 0
	    }
	    if { ! [info exists nonblocking(status,$count)] } {
		::tclu::ERROR "not a valid id, $count, for nonblocking exec ..."
		return 0
	    }	    
	    if { [info exists nonblocking(done,$count)] } {
		::tclu::ERROR "\"::tclu::nonblocking exec $count\" is locked, likely it has already been executed priviously"
		return 0
	    }
	    set args  [lrange $args 1 end]

	    set nonblocking(done,$count) 0
	    set nonblocking(fID,$count)  [open [concat | $args] r]
	    
	    fconfigure $nonblocking(fID,$count) -blocking 0
	    fileevent  $nonblocking(fID,$count) readable [list ::tclu::_nonblockingEvent $count]
	    
	    #tkwait variable ::tclu::nonblocking(done,$count)    
	    vwait ::tclu::nonblocking(done,$count)    
	    return $nonblocking(status,$count)
	}

	stdout {
	    set count [lindex $args 0]	    
	    set nonblocking(stdoutCmd,$count) [lindex $args 1]
	}
	
	save {
	    set count    [lindex $args 0]
	    set saveFile [lindex $args 1]
	    writeFile $saveFile $nonblocking(output,$count)
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
		if { ! [info exists nonblocking(done,$id)] } {
		    ::tclu::ERROR "nonblocking exec id $id does not exists"
		}
		catch {close $nonblocking(fID,$count)}
	    }
	}

	default {
	    ::tclu::usage "::tclu::nonblocking open|exec|unset|kill ?args?"
	}
    }
}

proc ::tclu::_nonblockingEvent {count} {
    variable nonblocking

    if { ! [eof $nonblocking(fID,$count)] } {

	set txt [gets $nonblocking(fID,$count)]
	append nonblocking(output,$count) $txt\n
	
	if { [info exists nonblocking(stdoutCmd,$count) ] } {
	    # form of stdout cmd is: 
	    #   either:  cmd id text
	    #   or:     {cmd arg1 arg2 ...} id text
	    # so we should use eval ...
	    eval $nonblocking(stdoutCmd,$count) $count [list $txt\n]
	}

    } else {
	if { [catch {close $nonblocking(fID,$count)}] } {
	    ::tclu::errorDialog "error while closing nonblocking execution \#.$count"
	    set nonblocking(status,$count) 0
	}

	if { [info exists nonblocking(stdoutCmd,$count) ] } {
	    #signal end of run to the stdout ...
	    eval $nonblocking(stdoutCmd,$count) $count [list {
		*** END-OF-RUN ***
	    }]
	}

	set nonblocking(done,$count) 1

    }
}



# ------------------------------------------------------------------------
#****f* ::tclu/::tclu::expandArgs
#  NAME
#    ::tclu::expandArgs -- expands the arguments
#  USAGE
#    tclu::expandArgs code
#  DESCRIPTION
#    This proc expands the argumants, for example:
#
#    The args may be in expanded or list form:
#        expanded-form:   -option value -option value
#        list-form:     { -option value -option value }
#
#    The the call \"tclu::expandArgs $args\" will return for both modes, 
#    the expanded from, i.e.: -option value -option value
#******
# ------------------------------------------------------------------------

proc ::tclu::expandArgs {code} {
    if { [llength $code] == 1 } {
        # list-mode
        return [lindex $code 0]
    } else {
        # expanded mode
        return $code
    }
}


# ------------------------------------------------------------------------
#****f* ::tclu/::tclu::getArgs
#  NAME
#    ::tclu::getArgs -- get requested option-value pairs from arguments
#  USAGE
#    tclu::getArgs optionList args
#  ARGUMENTS
#  * optionList -- list of selected options to return the option-value pairs
#  * args -- list of all option-value pairs
#  RETURN VALUE
#    The requested option-value pairs that were found in the arguments.
#******
# ------------------------------------------------------------------------

proc ::tclu::getArgs {optList args} {
    
    set args   [expandArgs $args]    
    set result ""
    
    foreach {opt value} $args {
        
        # is the option acceptable
        if { [lsearch -exact $optList $opt] > -1} {
            lappend result $opt $value
        }
    }
    
    return $result
}


# ------------------------------------------------------------------------
#****f* ::tclu/::tclu::extractArgs
# USAGE
#   tclu::extractArgs optionList argsVar
# ARGUMENTS
# * optionList -- list of options to extract
# * argsVar    -- name of the variable holding the arguments (option-value pairs)
# DESCRIPTION
#   Extract the requested options from arguments and remove them from
#   argsVar.
# RETURN VALUE
#   The requested option-value pairs that were found in the arguments.
#******
# ------------------------------------------------------------------------

proc ::tclu::extractArgs {optList argsVar} {
    upvar $argsVar args

    set result  ""
    set newArgs ""

    foreach {opt value} $args {        
        if { [lsearch -exact $optList $opt] > -1} {
            lappend result $opt $value
        } else {
	    lappend newArgs $opt
	    
	    if { $value != {} } {
		lappend newArgs $value
	    }
	}
    }
    set args $newArgs
    return $result
}


#****f* ::tclu/::tclu::modePatternArgs
# NAME
# modePatternArgs -- parse the args for mode & pattern
#
# USAGE
# tclu::modePatternArgs $modeVarname $ppaternVarname {$mode $pattern}
#
# DESCRIPTION
# This proc facilitates the use of $args as "mode pattern" in the
# style of "array names arrayName ?mode? ?pattern?".  If some proc
# accept $args in this style, then this routine helps. On exit the
# modeVar variable holds the mode, while patternVar the pattern.
#
# RETURN VALUE
# None.
#******
proc ::tclu::modePatternArgs {modeVar patternVar args} {
    upvar $modeVar mode
    upvar $patternVar pattern

    # args == {mode -glob} {namePattern *} 
    set mode        -glob
    set namePattern *

    set args [expandArgs $args]
    
    if { [llength $args] == 1 } {
	set pattern [lindex $args 0]
    } elseif { [llength $args] == 2 } {
	set mode    [lindex $args 0]
	set pattern [lindex $args 1]
    } elseif { [llength $args] > 2 } {
	error "wrong usage, should be get*** ?mode? ?pattern?"
    }
}

#****f* ::tclu/::tclu::mustBeEvenNumOfArgs
# NAME
#   mustBeEvenNumOfArgs -- trigger an error if number of arguments is not even
# USAGE
#   tclu::mustBeEvenNumOfArgs arguments
#******
proc ::tclu::mustBeEvenNumOfArgs {arg} {
    if { [llength $arg] % 2 } {
	uplevel {
	    error "number of arguments must be even, but the actual number is odd"
	}
    }
}


#****f* ::tclu/::tclu::mustBeOddNumOfArgs
# NAME
#   mustBeOddNumOfArgs -- trigger an error if number of arguments is not odd
# USAGE
#   tclu::mustBeOddNumOfArgs arguments
#******
proc ::tclu::mustBeOddNumOfArgs {arg} {
    if { ! ([llength $arg] % 2) } {	
	uplevel {
	    error "number of arguments must be odd, but the actual number is even"
	}
    }
}
