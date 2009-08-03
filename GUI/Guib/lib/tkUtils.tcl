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
# $Id: tkUtils.tcl,v 1.10 2009-08-03 14:22:43 kokalj Exp $ 
#

#------------------------------------------------------------------------
#****h* TclTkLib/::tku
#  NAME
#    TKU == Tone Kokalj's Tk Utilities
#    
#  COPYRIGHT
#    2002--2004 (c) by Tone Kokalj
#
#  AUTHOR
#    Tone Kokalj
#
#  CREATION DATE
#    Starting on Fri Dec 14 08:32:59 CET 2001
#
#  NOTES
#    Here are the basic Tk utility functions
#******
#------------------------------------------------------------------------

package provide tku 0.9
package require Tk

namespace eval ::tku {
    variable cursor 
    variable widgetCounter 0
    variable toplevelCounter 0
    variable getAllDescendants_list

    set cursor(default) [. cget -cursor]
    set cursor(watch)   watch
    
    namespace export widgetName
    namespace export toplevelName
    namespace export toplevelExists
    namespace export disableAll
    namespace export enableAll
    namespace export getAllDescendants
    namespace export centerWindow
    namespace export createFont
    namespace export setCursor
    namespace export resetCursor
    namespace export watchExec
    namespace export getOpenFile
    namespace export getSaveFile
    namespace export exitApp
}


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
#  *  parent   -- widget name of the parent widget
#  *  lastName -- a prefix for the lastname of the widget (i.e. .a.a.lastname)
#
#  RETURN VALUE
#    A unique widget name.
#
#  EXAMPLE
#    ::tku::widgetName $parent 
#********
#------------------------------------------------------------------------

proc ::tku::widgetName {{wid {}} {lastName {}}} {
    variable widgetCounter
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
#  * prefix   -- prefix for the widget name (for example ".tplw")
#
#  RETURN VALUE
#    A unique toplevel-widget name.
#
#  EXAMPLE
#    ::tku::toplevelName $parent 
#********
#------------------------------------------------------------------------

proc ::tku::toplevelName {{prefix .tplw}} {
    variable toplevelCounter

    set prefix .[string trim $prefix "."]

    while {1} { 
	if { [winfo exists ${prefix}$toplevelCounter] } {
	    incr toplevelCounter
	} else {
	    return ${prefix}$toplevelCounter
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
#    If it does not exists it just returns the pathName, otherwise the
#    "return -code return" is issued, to return from the caller proc.
#
#  ARGUMENTS
#  * pathName -- widget path-name
#
#  RETURN VALUE
#    If pathName widget does not exists it returns pathName, otherwise
#    returns the "-code return".
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



#****f* ::tku/::tku::disableAll
# SYNOPSIS
proc ::tku::disableAll {wlist} {
    # PURPOSE    
    #   Disables all wlist widgets and its children recursively.    
    # ARGUMENTS    
    # * wlist -- list of widgets to disable 
    # SOURCE

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

#******


#****f* ::tku/::tku::enableAll
# SYNOPSIS
proc ::tku::enableAll {wlist} {   
    # PURPOSE
    #  Enable all widgets and its children recursively. 
    # ARGUMENTS    
    # * wlist -- list of widgets to enable 
    # SOURCE

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
#******


#****f* ::tku/::tku::getAllDescendants
# SYNOPSIS
proc ::tku::getAllDescendants {w} {
    # PURPOSE    
    #   Return a list of all descendants of a given window.    
    # ARGUMENTS    
    # * w -- window for which to return all its descendants
    # SOURCE

    variable getAllDescendants_list
    
    if { [info exists getAllDescendants_list] } {
        set getAllDescendants_list ""
    }
    
    return [getAllDescendants_ $w]
}
proc getAllDescendants_ {wlist} {
    global getAllDescendants_list
    
    foreach w $wlist {  
        if { ![winfo exists $w] } continue
        
        lappend getAllDescendants_list $w
        
        set children [winfo children $w]
        
        if { $children != "" } {
            foreach child $children {
                getAllDescendants_ $child             
            }
        }
    }
    return $getAllDescendants_list
}
#******



#****f* ::tku/::tku::mouseWheelScrollCanvas
# SYNOPSIS
proc ::tku::mouseWheelScrollCanvas {can} {
    # PURPOSE 
    #   Arrange all necessary that a given canvas will be scrolled by
    #   mouse-wheel, even when the canvas possesses many descendant
    #   widgets (i.e. binds all descendants of a given canvas widget
    #   to scroll the canvas by mouse-wheel.
    # ARGUMENTS    
    # * can -- path of the canvas widget
    # SOURCE
    
    set scrollCmd_B4    [list $can yview scroll -2 units]
    set scrollCmd_B5    [list $can yview scroll +2 units]
    # TODO: please tune the %D on windows
    set scrollCmd_Wheel [list $can yview scroll %D units]
    
    set scrWin [getAllDescendants $can]
    
    global tcl_platform
    foreach wid $scrWin {
	if { [winfo exists $wid] } {
	    if { $tcl_platform(platform) == "unix" } {
		bind $wid <Button-4> $scrollCmd_B4
		bind $wid <Button-5> $scrollCmd_B5
	    } else {
		# TODO: please tune the %D on windows
		bind $wid <MouseWheel>  $scrollCmd_Wheel
	    }    
	}
    }
}
#******



#****f* ::tku/::tku::mouseWheelScrollDeleteBindings
# SYNOPSIS
proc ::tku::mouseWheelScrollDeleteBindings {can} {   
    # PURPOSE 
    #   The opposite of the ::tku::mouseWheelScrollCanvas command: it
    #   destroys all the mouse-wheel binding assciated with the
    #   scrolling of the given canvas, so that mouse-wheel does not
    #   have any effect more.
    # ARGUMENTS    
    # * can -- path of the canvas widget
    # SOURCE
    
    set scrWin [::tku::getAllDescendants $can]	
    global tcl_platform
    foreach wid $scrWin {
	if { [winfo exists $wid] } {
	    if { $tcl_platform(platform) == "unix" } {
		bind $wid <Button-4> ""
		bind $wid <Button-5> ""
	    } else {
		# TODO: please tune the %D on windows
		bind $wid <MouseWheel>  ""
	    }    
	}
    }
}    



#****f* ::tku/::tku::centerWindow
# SYNOPSIS
proc ::tku::centerWindow {thisWin {otherWid {}}} {
    # PURPOSE
    #   Centers the toplevel with respect to another widget or the
    #   screen as a whole.    
    # ARGUMENTS
    # * thisWin -- name of widget to center
    # * otherWid -- name of widhet the thisWin will be centered onto
    # SOURCE

    update idletasks

    # if otherWid is {}, then center width respect to the root window
 
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
        set reqX [expr {[winfo rootx $otherWid]+($otherWidW-($otherWidW/2))-($w/2)}]
        set reqY [expr {[winfo rooty $otherWid]+($otherWidH-($otherWidH/2))-($h/2)}]

        
        # Adjust for errors - if too long or too tall
        
        if { $reqX+$w+$wfudge > $rw } { set reqX [expr {$rw-$w-$wfudge}] }
        if { $reqX < $wfudge } { set reqX $wfudge }
        if { $reqY+$h+$hfudge > $rh } { set reqY [expr {$rh-$h-$hfudge}] }
        if { $reqY < $hfudge } { set reqY $hfudge }
    } 
    
    wm geometry $thisWin +$reqX+$reqY
}
#******


#****f* ::tku/::tku::createFont
# SYNOPSIS
proc ::tku::createFont {args} {
    # PURPOSE
    #  Create new Tk Font with requested attributes (args=="option value"
    #  pairs), and return its name
    # ARGUMENTS
    # * args -- arguments passed to font configure (i.e. option value
    #   pairs specifying font attributres, see font Tk-command) 
    # 
    # SOURCE

    set fontName [font create]
    eval font configure $fontName $args
    return $fontName
}
#******



#****f* ::tku/::tku::setCursor
# SYNOPSIS
proc ::tku::setCursor {cursorName {window .}} {

    # PURPOSE
    #   Set a new cursor, one among those defined by cursor variables.
    #   So far supported names are: 
    #     * watch
    #     * default
    # ARGUMENTS
    # * cursorName -- name for the new cursor
    # * window -- (optional) name of window for which to change cursor
    #   (propagates to its children too) 
    # SOURCE

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
#******



#****f* ::tku/::tku::resetCursor
# SYNOPSIS
proc ::tku::resetCursor {{window .}} {
    # PURPOSE
    #   Sets a default cursor for specified window.    
    # ARGUMENTS    
    # * window -- (optional) name of window for which to set the cursor to default
    #   (propagates to its children too) 
    # SOURCE

    variable cursor

    foreach t [winfo children $window] {
	if { [info commands $t] != "" } {
	    $t config -cursor $cursor(default)
	}
    }
    $window config -cursor $cursor(default)
    update
}
#******


#------------------------------------------------------------------------
#****f* ::tku/::tku::watchExec
#  NAME
#    ::tku::watchExec -- execute a script and display a watch cursor
#  USAGE
#    ::tku::watchExec script
#
#  DESCRIPTION
#    This proc executes a script and displays a watch cursor during exec. The
#    content of the script is evaluated at the caller level (i.e. uplevel 1
#
#  ARGUMENTS
#  *  parent   -- widget name of the parent widget
#  *  lastName -- a prefix of the last-name of the widget (i.e. .a.a.lastname)
#
#  RETURN VALUE
#    Returns the last return-value of the last command in the script.
#
#  EXAMPLE
#    ::tku::watchExec { 
#        do_whatever $a $b 
#        do_this_and_that $c $d
#    }
#********
#------------------------------------------------------------------------

proc ::tku::watchExec {script} {
    ::tku::setCursor watch
    set result [uplevel 1 $script]
    ::tku::resetCursor
    return $result
}



#****f* ::tku/::tku::getOpenFile
# SYNOPSIS
proc ::tku::getOpenFile {file args} {
    # PURPOSE
    #   Get (return) the name of file for opening, in particular:
    #      if $file == "" --> query file
    #      else           --> check if file exists
    # ARGUMENTS
    # * file -- name of file to open (if empty-string, query the file to open)
    # * args -- arguments passed to tk_getOpenFile if file == {}
    # SOURCE

    if { $file == {} } {
	
	# query file ...
	
	set file [eval tk_getOpenFile $args]
	if { $file == "" } {
	    return -code return
	}
    } else {
	
	# check if file exists
	
	if { ![file exists $file] } {
	    # note: tku's errorDialog should not use return
	    ::tclu::warningDialog "file \"$file\" does not exists !!!"
	    return -code return
	}
    }

    return $file
}
#******



#****f* ::tku/::tku::getSaveFile
# SYNOPSIS
proc ::tku::getSaveFile {file args} {
    # PURPOSE
    #   Get (return) the name of file for saving, in particular:
    #      if $file == "" --> query file
    #      else           --> check if dirname of file exists and is writeable
    # ARGUMENTS
    # * file -- name of file to open (if empty-string, query the file to open)
    # * args -- arguments passed to tk_getSaveFile if file == {}
    # SOURCE

    if { $file == {} } {
		
	# query file ...
	
	set file [eval tk_getSaveFile $args]
	if { $file == "" } {
	    return -code return
	}
    } else {
	
	# check if dirname of file exists
	
	set dirname [file dirname $file]
	if { ! [file writable $dirname] } {
	    ::tclu::errorDialog "can't create file \"$file\". Permission denied !!!"
	    return -code return ""
	}
	if { ! [file isdirectory $dirname] } {
	    # note: tku's errorDialog should not use return
	    tclu::errorDialog "can't create \"$file\". Directory $dirname does not exists. !!!"
	    return -code return ""
	}
    }
    
    return $file
}
#******



#****f* ::tku/::tku::exitApp
# SYNOPSIS
proc ::tku::exitApp {{cmdAtExit {}}} {
    # PURPOSE
    #   Exit from application if user clicks "Yes" button. 
    # ARGUMENTS    
    # * cmdAtExit -- command to execute before exiting 
    # SOURCE

    set button [tk_messageBox -message "Really quit?" \
                    -type yesno -icon question]
    if { $button == "yes" } {
        # do everything requested before exiting ... 
	if { $cmdAtExit != "" } {
	    eval $cmdAtExit
	}

        # delete all temporary files ...
        ::tclu::tempFile delete all
        exit
    }
}
#******
