#
# $RCSfile: tkUtils.tcl,v $ --
#
#      This file contains the Tone Kokalj's Tk utilities functions.
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
# $Id: tkUtils.tcl,v 1.1 2004-02-18 11:29:17 kokalj Exp $ 
#

#------------------------------------------------------------------------
#****h* TkLib/::tku
#  NAME
#    TKU == Tone Kokalj's Tk Utilities
#                         ^^ ^ == tku
#  COPYRIGHT
#    2002--2004 (c) by Tone Kokalj
#
#  AUTHOR
#    Tone Kokalj
#
#  CREATION DATE
#    Starting on Fri Dec 14 08:32:59 CET 2001
#
#  MODIFICATION HISTORY
#    Writing a few functions (based on XCRYSDEN's Tcl functions)
#
#  NOTES
#    Here are the basic Tcl-only utility functions
#******
#------------------------------------------------------------------------

package provide tku 0.1

namespace eval ::tku:: {
    variable cursor 

    set cursor(default) [. cget -cursor]
    set cursor(watch)   watch

    #::tku::_init
    
    namespace export widgetName
    namespace export toplevelName
    namespace export toplevelExists
    namespace export disableAll
    namespace export enableAll
    namespace export errorDialog
    namespace export centerWindow
    namespace export createFont
    namespace export setCursor
    namespace export resetCursor
}

#proc ::tku::_init {} {
#    variable cursor
#    set cursor(default) [. cget -cursor]
#    set cursor(watch)   watch
#}
   


#------------------------------------------------------------------------
#****f* ::tku/::tku::widgetName
#  NAME
#    ::tku::widgetName -- generates new widget name
#  USAGE
#    ::tku::widgetName ?parent? ?lastName?
#
#  DESCRIPTION
#    This proc is used for generating a unique widget name.
#
#  ARGUMENTS
#    parent   -- widget name of the parent widget
#    lastName -- a prefix of the last-name of the widget (i.e. .a.a.lastname)
#
#  RETURN VALUE
#    A unique widget name.
#
#  EXAMPLE
#    ::tku::widgetName $parent 
#********
#------------------------------------------------------------------------

proc ::tku::widgetName {{wid {}} {lastName {}}} {
    set widgetCounter 0
    while {1} { 
	if { [winfo exists ${wid}.${lastName}$widgetCounter] } {
	    incr widgetCounter
	} else {
	    return ${wid}.${lastName}$widgetCounter
	}
    }
}

#------------------------------------------------------------------------
#****f* ::tku/::tku::toplevelName
#  NAME
#    ::tku::toplevelName -- generates new toplevel-widget name
#  USAGE
#    ::tku::toplevelName ?prefix?
#
#  DESCRIPTION
#    This proc is used for generating a unique toplevel-widget name.
#
#  ARGUMENTS
#    prefix   -- prefix for the widget name (for example ".tplw")
#
#  RETURN VALUE
#    A unique toplevel-widget name.
#
#  EXAMPLE
#    ::tku::toplevelName $parent 
#********
#------------------------------------------------------------------------

proc ::tku::toplevelName {{prefix .tplw}} {
    set prefix .[string trim $prefix "."]
    set widgetCounter 0
    while {1} { 
	if { [winfo exists ${prefix}$widgetCounter] } {
	    incr widgetCounter
	} else {
	    return ${prefix}$widgetCounter
	}
    }
}


#------------------------------------------------------------------------
#****f* ::tku/::tku::toplevelExists
#  NAME
#    ::tku::toplevelExists -- checks if toplevel already exists
#  USAGE
#    ::tku::toplevelExists pathName
#
#  DESCRIPTION
#    This proc is used for checking if widget pathName already exists.
# If it does not existsit just returns the pathName, otherwise the
# "return -code return" is issued, to return from the caller proc.

#
#  ARGUMENTS
#    pathName -- widget path-name
#
#  RETURN VALUE
#    If pathName widget does not exists it returns pathName, otherwise
#    returns the -code return.
#
#  EXAMPLE
#    ::tku::toplevelExists $widget
#********
#------------------------------------------------------------------------

proc ::tku::toplevelExists {pathName} {
    if { [winfo exists $pathName] } {
	return -code return
    } else {
	return $pathName
    }
}


# ------------------------------------------------------------------------
#  disables all wlist widgets and its children recursively
# ------------------------------------------------------------------------
proc ::tku::disableAll {wlist} {
    foreach w $wlist {
	if { ! [winfo exists $w] } { 
	    continue 
	}
	set children [winfo children $w]
	if { $children != "" } {
	    foreach w $children {
		disableAll $w
		catch {$w configure -state disabled}
	    }
	}
    }	
}

# ------------------------------------------------------------------------
#  enable all widgets and its children recursively
# ------------------------------------------------------------------------
proc ::tku::enableAll {wlist} {   
    foreach w $wlist {
	if { ! [winfo exists $w] } {
	    continue
	}
	set children [winfo children $w]
	if { $children != "" } {
	    foreach w $children {
		enableAll $w
		catch {$w configure -state normal}
	    }
	}
    }
}


# ------------------------------------------------------------------------
#  print error message and returns from the caller proc
# ------------------------------------------------------------------------
proc ::tku::errorDialog {text} {
    tk_messageBox -title ERROR -message "ERROR: $text" -type ok -icon error
    return -code return ""
}

# ------------------------------------------------------------------------
#  Centers the toplevel with respect to another widget or the screen
#  as a whole. 
# ------------------------------------------------------------------------
proc ::tku::centerWindow {thisWin {otherWid {}}} {
    update idletasks

    # if otherWid is {} center width respect to the root window
 
    set w  [winfo reqwidth $thisWin]
    set h  [winfo reqheight $thisWin]
    # root window height/width
    set rh [winfo screenheight $thisWin]     
    set rw [winfo screenwidth $thisWin]
 
    if { $otherWid == "" } {
        set reqX [expr {($rw-$w)/2}]
        set reqY [expr {($rh-$h)/2}]
    } else {
        set wfudge 5      ;# wm width fudge factor
        set hfudge 20     ;# wm height fudge factor
        set otherWidW [winfo width $otherWid]
        set otherWidH [winfo height $otherWid]
        set reqX [expr [winfo rootx \
			    $otherWid]+($otherWidW-($otherWidW/2))-($w/2)]
        set reqY [expr [winfo rooty \
			    $otherWid]+($otherWidH-($otherWidH/2))-($h/2)]

        #
        # Adjust for errors - if too long or too tall
        #
        if { [expr $reqX+$w+$wfudge] > $rw } { set reqX [expr $rw-$w-$wfudge] }
        if { $reqX < $wfudge } { set reqX $wfudge }
        if { [expr $reqY+$h+$hfudge] > $rh } { set reqY [expr $rh-$h-$hfudge] }
        if { $reqY < $hfudge } { set reqY $hfudge }
    } 
    
    wm geometry $thisWin +$reqX+$reqY
}

# ------------------------------------------------------------------------
#  Create new Tk Font with requested attributes (args=="option value"
#  pairs), and return its name
# ------------------------------------------------------------------------
proc ::tku::createFont {args} {
    set fontName [font create]
    eval font configure $fontName $args
    return $fontName
}


proc ::tku::setCursor {cursorName {window .}} {
    variable cursor

    switch -exact -- $cursorName {
	watch   { set _cursor $cursor(watch)   }
	default { set _cursor $cursor(default) }
    }

    foreach t [winfo children $window] {
	if { [info commands $t] != "" } {
	    $t config -cursor $_cursor
	}
    }
    $window config -cursor $_cursor
    update
}

proc ::tku::resetCursor {{window .}} {
    variable cursor

    foreach t [winfo children $window] {
	if { [info commands $t] != "" } {
	    $t config -cursor $cursor(default)
	}
    }
    $window config -cursor $cursor(default)
    update
}


    
